using System;
using System.Collections;
using System.Collections.Generic;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

// <Custom using>
using System.Drawing;
using System.Threading;
using System.Threading.Tasks;
using System.Collections.Concurrent;
// </Custom using>



/// <summary>
/// This class will be instantiated on demand by the Script component.
/// </summary>
public class Script_Instance : GH_ScriptInstance
{
    #region Utility functions
    /// <summary>Print a String to the [Out] Parameter of the Script component.</summary>
    /// <param name="text">String to print.</param>
    private void Print(string text) { __out.Add(text); }
    /// <summary>Print a formatted String to the [Out] Parameter of the Script component.</summary>
    /// <param name="format">String format.</param>
    /// <param name="args">Formatting parameters.</param>
    private void Print(string format, params object[] args) { __out.Add(string.Format(format, args)); }
    /// <summary>Print useful information about an object instance to the [Out] Parameter of the Script component. </summary>
    /// <param name="obj">Object instance to parse.</param>
    private void Reflect(object obj) { __out.Add(GH_ScriptComponentUtilities.ReflectType_CS(obj)); }
    /// <summary>Print the signatures of all the overloads of a specific method to the [Out] Parameter of the Script component. </summary>
    /// <param name="obj">Object instance to parse.</param>
    private void Reflect(object obj, string method_name) { __out.Add(GH_ScriptComponentUtilities.ReflectType_CS(obj, method_name)); }
    #endregion

    #region Members
    /// <summary>Gets the current Rhino document.</summary>
    private RhinoDoc RhinoDocument;
    /// <summary>Gets the Grasshopper document that owns this script.</summary>
    private GH_Document GrasshopperDocument;
    /// <summary>Gets the Grasshopper script component that owns this script.</summary>
    private IGH_Component Component;
    /// <summary>
    /// Gets the current iteration count. The first call to RunScript() is associated with Iteration==0.
    /// Any subsequent call within the same solution will increment the Iteration count.
    /// </summary>
    private int Iteration;
    #endregion

    /// <summary>
    /// This procedure contains the user code. Input parameters are provided as regular arguments, 
    /// Output parameters as ref arguments. You don't have to assign output parameters, 
    /// they will have a default value.
    /// </summary>
    private void RunScript(bool reset, bool go, bool debug, int method, Mesh M, List<Point3d> P, double sR, double sI, double dR, double eR, double sensDist, double sensAng, ref object points, ref object vectors, ref object outMCol, ref object neigh, ref object neighCol, ref object sensors)
    {
        // <Custom code> 
        // ............................................................................
        // Stigmergy system
        /*
         notes on code refactoring

        build 2 main classes:
        . Particle System (and a class Particle)
        . Environment (Mesh with scalar and optional vector field)

        RunScript function structure at the end:

        . check inputs & eventually bypass execution
        . initialise Environment
        . initialise Paeticle System
        . Update live variables
        . Update Particle System
        . Update Environment
        . Extract Output Geometries and Data

         */

        // return on null input
        if (M == null || P == null) return;


        // initialize on first execution or mesh change
        if (PS == null || M.Vertices.Count != scalarField.Length)
        {
            // initialize particle system
            PS = new ParticleSystem(P, M);

            // populate scalar field from mesh and array to restore it
            scalarField = PopulateScalarField(M);
            initVal = scalarField;

            // populate point cloud from mesh (speeds up neighbour search)
            MeshPoints = PopulatePointCloud(M);

            // build neighbours indexes map & diffusion weights
            NeighbourMap = BuildNeighboursMap(M);
            DiffusionWeights = CalculateDiffusionWeights();
        }

        // restore initial values on reset
        if (reset)
        {
            // initialize particle system
            PS = new ParticleSystem(P, M);

            // restore scalar field values
            scalarField = RestoreScalarField(initVal);

        }

        if (go)
        {
            // update runtime variables
            PS.seekRadius = sR;
            PS.seekIntensity = sI;
            PS.sensDist = sensDist; // old value 3.0
            PS.sensAng = sensAng * 2;

            // update simulation
            switch (method)
            {
                case 0:
                    PS.UpdateRTree(M);
                    break;
                case 1:
                    PS.UpdateJones();
                    break;
            }

            // diffusion
            Diffusion((float)dR);

            // evaporate scalar field
            EvaporateField(scalarField, eR);

            // update component
            Component.ExpireSolution(true);
        }

        // . . . . . . .  extract geometries

        // particles positions and velocites
        PS.GetPointsVectors(out pts, out vecs);
        points = pts;
        vectors = vecs;
        // colored mesh
        outMCol = GetColoredMesh(M, scalarField);
        // debug mode
        if (debug)
        {
            switch (method)
            {
                case 0:
                    neigh = PS.GetNeighPts();
                    neighCol = PS.GetNeighBrightness();
                    break;
                case 1:
                    sensors = PS.SensorsOut();
                    break;
            }
        }

        // ............................................................................
        // </Custom code> 
    }

    // <Custom additional code> 
    // ............................................................................Global Variables
    public ParticleSystem PS;
    public MeshEnvironment ME;
    public GH_Point[] pts;
    public GH_Vector[] vecs;
    public static float[] scalarField;
    public static PointCloud MeshPoints;
    public static int[][] NeighbourMap;
    public static float[] DiffusionWeights;
    public float[] initVal;
    public static DataTree<Vector3d> ParticleSensors;

    // ............................................................................Classes

    // .......................................................... Particle System ........

    public class ParticleSystem
    {
        // . . . . . . . . . . . . . . . . . . . . . . fields
        public List<Particle> Particles;
        public RTree MeshRTree;
        public double seekRadius;
        public double seekIntensity;
        public double sensDist;
        public double sensAng;

        // . . . . . . . . . . . . . . . . . . . . . . constructor
        public ParticleSystem(List<Point3d> positions, Mesh M)
        {
            Particles = new List<Particle>();
            for (int i = 0; i < positions.Count; i++)
                Particles.Add(new Particle(this, positions[i], RandomVectorOnMesh(M, positions[i], i) * 1.5));

            MeshRTree = PopulateMeshRTree(M);
        }

        // . . . . . . . . . . . . . . . . . . . . . . methods
        /// <summary>
        /// Update with Jeff Jones 3-sensors method
        /// </summary>
        public void UpdateJones()
        {
            /*
             How was it parallelized?
             Several particles might write their pheromone contribution to the same mesh point,
             to avoid simultaneous write to te same location (surest way to get an error,
             even using so called "thread-safe" collections - especially if they are made
             of structures and not classes):
             . each particle saves its own closest point index as an internal field
             . the resulting contribution is written in a separate non-parallel loop
               that updates the scalar field.
             */
            Parallel.ForEach(Particles, p =>

            {
                // find mesh closest point index
                p.currentCPindex = MeshPoints.ClosestPoint(p.pos);
                // Seek brightest point direction
                p.SeekAngVis(seekIntensity);

            });

            // write to scalar field
            foreach (Particle p in Particles)
                scalarField[p.currentCPindex] = (float)Math.Min(1.0, scalarField[p.currentCPindex] + p.PheroStrength);

            // update particles
            if (Particles.Count > 50)
            {
                Parallel.ForEach(Particles, p =>
                {
                    p.Update();
                });
            }
            else
                foreach (Particle p in Particles)
                    p.Update();

        }

        /// <summary>
        /// Update with RTree method (samples points within a sphere)
        /// </summary>
        /// <param name="M"></param>
        public void UpdateRTree(Mesh M)
        {
            foreach (Particle p in Particles)
            {
                // clear particle neighbours list
                p.neighbours.Clear();
                p.neighFieldPts.Clear();

                // Eventhandler function for RTree search
                EventHandler<RTreeEventArgs> rTreeCallback = (object sender, RTreeEventArgs args) =>
                {
                    p.neighbours.Add(M.Vertices[args.Id]);
                    p.neighFieldPts.Add(new FieldPt(args.Id, scalarField[args.Id]));
                };

                MeshRTree.Search(new Sphere(M.ClosestPoint(p.pos + p.vel * sensDist), seekRadius), rTreeCallback);

                // Seek brightest point direction
                p.SeekColor(seekIntensity);

                //// update scalar field
                int cpInd = MeshPoints.ClosestPoint(p.pos);
                scalarField[cpInd] = (float)Math.Min(1.0, scalarField[cpInd] + p.PheroStrength);

            }

            if (Particles.Count > 50)
            {
                Parallel.ForEach(Particles, p =>
               {
                   p.Update();
               });
            }
            else
                foreach (Particle p in Particles)
                    p.Update();

        }

        public void GetPointsVectors(out GH_Point[] pts, out GH_Vector[] vecs)
        {
            pts = new GH_Point[Particles.Count];
            vecs = new GH_Vector[Particles.Count];

            for (int i = 0; i < Particles.Count; i++)
            {
                pts[i] = new GH_Point(Particles[i].pos);
                vecs[i] = new GH_Vector(Particles[i].vel);
            }
        }

        public DataTree<GH_Point> GetNeighPts()
        {
            DataTree<GH_Point> ptsOut = new DataTree<GH_Point>();

            for (int i = 0; i < Particles.Count; i++)
                ptsOut.EnsurePath(new GH_Path(i));

            // parallelize for more than 50 particles
            if (Particles.Count > 50)
            {

                var x = Partitioner.Create(0, Particles.Count);

                Parallel.ForEach(x, (range, loopstate) =>
                {
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        for (int j = 0; j < Particles[i].neighbours.Count; j++)
                            ptsOut.Add(new GH_Point(Particles[i].neighbours[j]), new GH_Path(i));
                    }
                });
            }
            else
            {
                for (int i = 0; i < Particles.Count; i++)
                {
                    for (int j = 0; j < Particles[i].neighbours.Count; j++)
                        ptsOut.Add(new GH_Point(Particles[i].neighbours[j]), new GH_Path(i));
                }
            }

            return ptsOut;
        }

        public DataTree<GH_Number> GetNeighBrightness()
        {

            DataTree<GH_Number> briOut = new DataTree<GH_Number>();

            for (int i = 0; i < Particles.Count; i++)
                briOut.EnsurePath(new GH_Path(i));

            for (int i = 0; i < Particles.Count; i++)
            {
                for (int j = 0; j < Particles[i].neighFieldPts.Count; j++)
                    briOut.Add(new GH_Number(Particles[i].neighFieldPts[j].scalar), new GH_Path(i));
            }

            return briOut;
        }

        public DataTree<GH_Vector> SensorsOut()
        {
            DataTree<GH_Vector> sensOut = new DataTree<GH_Vector>();

            for (int i = 0; i < Particles.Count; i++)
            {
                GH_Path path = new GH_Path(i);
                sensOut.Add(new GH_Vector(Particles[i].sensorL), path);
                sensOut.Add(new GH_Vector(Particles[i].sensorC), path);
                sensOut.Add(new GH_Vector(Particles[i].sensorR), path);
            }
            return sensOut;
        }

    }

    // .......................................................... Particle ...............
    public class Particle
    {
        public Point3d pos;
        public Vector3d vel;
        public Vector3d acc;
        public Vector3d sensorC, sensorL, sensorR;
        public double visAng;
        public double MaxSpeed;
        public float PheroStrength;
        public int currentCPindex;
        public List<Point3d> neighbours;
        public List<FieldPt> neighFieldPts;
        public ParticleSystem pSystem;

        public Particle(ParticleSystem pSystem, Point3d pos, Vector3d vel)
        {
            this.pSystem = pSystem;
            this.pos = pos;
            this.vel = vel;
            acc = Vector3d.Zero;
            MaxSpeed = 1.5f;
            PheroStrength = 0.1f;
            visAng = 0.3 * Math.PI;
            neighbours = new List<Point3d>();
            neighFieldPts = new List<FieldPt>();
        }

        public void Update()
        {
            Move();
        }

        public void SeekAngVis(double seekIntensity)
        {
            acc = Vector3d.Zero;
            visAng = pSystem.sensAng;
            Vector3d rotAxis = MeshPoints[currentCPindex].Normal;
            sensorC = VectorAmp(vel, pSystem.sensDist);

            sensorL = new Vector3d(sensorC);
            sensorR = new Vector3d(sensorC);
            sensorL.Rotate(visAng * 0.5, rotAxis);
            sensorR.Rotate(visAng * -0.5, rotAxis);

            int pL, pR, pC;
            double briL, briR, briC;
            pC = MeshPoints.ClosestPoint(pos + sensorC);
            pL = MeshPoints.ClosestPoint(pos + sensorL);
            pR = MeshPoints.ClosestPoint(pos + sensorR);
            briC = scalarField[pC];
            briL = scalarField[pL];
            briR = scalarField[pR];
            Vector3d desired = Vector3d.Zero;

            //if (briL - briR < 0.01)
            //{
            //    Random rnd = new Random();
            //    futPos.Rotate(visAng * (rnd.NextDouble() - 0.5), rotAxis);
            //}
            //else
            //    futPos = briL > briR ? futPosL : futPosR;

            // find brightest sensor
            Vector3d direction = sensorC;
            if (briL > briC && briL > briR) direction = sensorL;
            else if (briR > briC) direction = sensorR;

            // find mesh closest point to futpos
            desired = MeshPoints[MeshPoints.ClosestPoint(pos + direction)].Location - pos;
            desired.Unitize();
            desired *= MaxSpeed;

            acc = (desired - vel) * seekIntensity;
        }

        public void SeekColor(double seekIntensity)
        {
            acc = vel;

            if (neighFieldPts.Count > 0)
            {
                acc = Vector3d.Zero;
                Vector3d desired = Vector3d.Zero;
                double bri;
                double MaxBri = -1.0;
                int MaxInd = -1;

                for (int i = 0; i < neighFieldPts.Count; i++)
                {
                    bri = neighFieldPts[i].scalar; // brightness is 0-1 and RGB 0-255
                    if (bri > MaxBri)
                    {
                        MaxBri = bri;
                        MaxInd = i;
                    }

                }
                desired += neighbours[MaxInd] - pos;

                desired.Unitize();
                desired *= MaxSpeed;

                acc = (desired - vel) * seekIntensity;
            }
        }

        //public void SeekMeshClosestPt(Mesh M, double seekIntensity)
        //{
        //    acc = Vector3d.Zero;
        //    MeshPoint mP = M.ClosestMeshPoint(pos + vel * pSystem.futPosMult, 100); // old value 1.5
        //    if (mP != null)
        //    {
        //        Point3d mCP = mP.Point;

        //        Vector3d desired = mCP - pos;
        //        desired.Unitize();
        //        desired *= MaxSpeed;

        //        acc = (desired - vel) * seekIntensity;
        //    }

        //}

        public void Move()
        {
            //vel = vel * 0.9 + acc * 0.1;
            vel = vel + acc;
            if (vel.Length > MaxSpeed)
            {
                vel.Unitize();
                vel *= MaxSpeed;
            }
            if (vel.Length < 1)
            {
                vel.Unitize();
            }
            pos += vel;
        }


    }

    // .......................................................... Mesh Environment .......

    public class MeshEnvironment
    {
        public Mesh M;
        public float diffusionRate;
        public float evaporationRate;
        public float[] scalarField;
        float[] initVal;
        public PointCloud MeshPoints;
        public RTree MeshRTree;
        int[][] NeighbourMap;
        float[] DiffusionWeights;
        

        // M must be a colored Mesh
        public MeshEnvironment(Mesh M, double diffusionRate, double evaporationRate)
        {
            // assign external variables
            this.M = M;
            this.diffusionRate = (float)diffusionRate;
            this.evaporationRate = (float)evaporationRate;
            // scalar field structures
            scalarField = PopulateScalarField(M);
            initVal = scalarField;
            NeighbourMap = BuildNeighboursMap(M);
            DiffusionWeights = CalculateDiffusionWeights();
            // point data structures
            MeshPoints = PopulatePointCloud(M);
            MeshRTree = PopulateMeshRTree(M);
        }

        public void Update()
        {
            Diffusion();
            EvaporateField();
        }

        public float[] PopulateScalarField(Mesh M)
        {
            float[] scalarField = new float[M.Vertices.Count];

            if (M.VertexColors.Count != 0)
            {

                var x = Partitioner.Create(0, M.VertexColors.Count);

                Parallel.ForEach(x, (range, loopstate) =>
                {
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        scalarField[i] = M.VertexColors[i].GetBrightness();

                    }
                });
            }
            return scalarField;
        }

        public void EvaporateField()
        {

            Parallel.For(0, scalarField.Length, i =>
            {
                scalarField[i] *= evaporationRate;
            });
        }

        public void RestoreScalarField()
        {
            Parallel.For(0, initVal.Length, i =>
            {
                scalarField[i] = initVal[i];
            });
        }

        public int[][] BuildNeighboursMap(Mesh M)
        {
            int[][] NeighbourMap = new int[M.Vertices.Count][];

            Parallel.For(0, M.Vertices.Count, i =>
            {
                NeighbourMap[i] = M.Vertices.GetConnectedVertices(i);
            });

            return NeighbourMap;

        }

        public float[] CalculateDiffusionWeights()
        {
            float[] Weights = new float[NeighbourMap.Length];

            Parallel.For(0, NeighbourMap.Length, i =>
            {
                Weights[i] = 1 / (float)NeighbourMap[i].Length;
            });

            return Weights;
        }

        public float[] Diffusion()
        {
            float[] newVal = new float[scalarField.Length];

            // calculate new scalar values
            Parallel.For(0, scalarField.Length, i =>
            {
                float neighVal = 0;

                for (int j = 0; j < NeighbourMap[i].Length; j++)
                    neighVal += scalarField[NeighbourMap[i][j]] * DiffusionWeights[i];
                newVal[i] = scalarField[i] * (1 - diffusionRate) + neighVal * diffusionRate;
            });

            // assign new scalar values
            Parallel.For(0, scalarField.Length, i =>
            {
                scalarField[i] = newVal[i];
            });

            return scalarField;
        }

        public RTree PopulateMeshRTree(Mesh M)
        {
            RTree rt = new RTree();

            for (int i = 0; i < M.Vertices.Count; i++)
            {
                rt.Insert((Point3d)M.Vertices[i], i);
            }

            return rt;
        }

        public PointCloud PopulatePointCloud(Mesh M)
        {
            PointCloud mP = new PointCloud();

            for (int i = 0; i < M.Vertices.Count; i++)
            {
                mP.Add(M.Vertices[i], M.Normals[i]);
            }

            return mP;
        }

        // output data

        public GH_Mesh GetColoredMesh()
        {

            Parallel.For(0, scalarField.Length, i =>
            {
                int c = (int)(scalarField[i] * 255);
                M.VertexColors[i] = Color.FromArgb(c, c, c);
            });

            return new GH_Mesh(M);
        }

        public GH_Number[] GetScalarField()
        {
            GH_Number[] outField = new GH_Number[scalarField.Length];

            Parallel.For(0, scalarField.Length, i =>
            {
                outField[i] = new GH_Number(scalarField[i]);
            });

            return outField;
        }
    }

    // .......................................................... Field Point ............
    public class FieldPt
    {
        public int id;
        public float scalar;
        public Vector3d vector;

        public FieldPt(int id, float scalar, Vector3d vector)
        {
            this.id = id;
            this.scalar = scalar;
            this.vector = vector;
        }

        public FieldPt(int id, float scalar) : this(id, scalar, Vector3d.Zero) { }

    }

    // ............................................................................Utilities

    #region Utilities

    public GH_Mesh GetColoredMesh(Mesh M, float[] scalarField)
    {

        Parallel.For(0, scalarField.Length, i =>
        {
            int c = (int)(scalarField[i] * 255);
            M.VertexColors[i] = Color.FromArgb(c, c, c);
        });

        return new GH_Mesh(M);
    }

    public GH_Number[] GetScalarField(float[] scalarField)
    {
        GH_Number[] outField = new GH_Number[scalarField.Length];

        Parallel.For(0, scalarField.Length, i =>
        {
            outField[i] = new GH_Number(scalarField[i]);
        });

        return outField;
    }

    public static void EvaporateField(float[] scalarField, double evapoRatio)
    {
        float evap = (float)evapoRatio;

        Parallel.For(0, scalarField.Length, i =>
        {
            scalarField[i] *= evap;
        });
    }

    public float[] PopulateScalarField(Mesh M)
    {
        float[] scalarField = new float[M.Vertices.Count];

        if (M.VertexColors.Count != 0)
        {

            var x = Partitioner.Create(0, M.VertexColors.Count);

            Parallel.ForEach(x, (range, loopstate) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    scalarField[i] = M.VertexColors[i].GetBrightness();

                }
            });
        }
        return scalarField;
    }

    public int[][] BuildNeighboursMap(Mesh M)
    {
        int[][] NeighbourMap = new int[M.Vertices.Count][];

        Parallel.For(0, M.Vertices.Count, i =>
        {
            NeighbourMap[i] = M.Vertices.GetConnectedVertices(i);
        });

        return NeighbourMap;

    }

    public float[] CalculateDiffusionWeights()
    {
        float[] Weights = new float[NeighbourMap.Length];

        Parallel.For(0, NeighbourMap.Length, i =>
        {
            Weights[i] = 1 / (float)NeighbourMap[i].Length;
        });

        return Weights;
    }

    public float[] Diffusion(float diffusionRate)
    {
        float[] newVal = new float[scalarField.Length];

        // calculate new scalar values
        Parallel.For(0, scalarField.Length, i =>
        {
            float neighVal = 0;

            for (int j = 0; j < NeighbourMap[i].Length; j++)
                neighVal += scalarField[NeighbourMap[i][j]] * DiffusionWeights[i];
            newVal[i] = scalarField[i] * (1 - diffusionRate) + neighVal * diffusionRate;
        });

        // assign new scalar values
        Parallel.For(0, scalarField.Length, i =>
        {
            scalarField[i] = newVal[i];
        });

        return scalarField;
    }

    public float[] RestoreScalarField(float[] initVal)
    {
        float[] scalarField = new float[initVal.Length];

        Parallel.For(0, initVal.Length, i =>
        {
            scalarField[i] = initVal[i];
        });

        return scalarField;
    }

    public static RTree PopulateMeshRTree(Mesh M)
    {
        RTree rt = new RTree();

        for (int i = 0; i < M.Vertices.Count; i++)
        {
            rt.Insert((Point3d)M.Vertices[i], i);
        }

        return rt;
    }

    public PointCloud PopulatePointCloud(Mesh M)
    {
        PointCloud mP = new PointCloud();

        for (int i = 0; i < M.Vertices.Count; i++)
        {
            mP.Add(M.Vertices[i], M.Normals[i]);
        }

        return mP;
    }

    public static Vector3d VectorAmp(Vector3d v, double Amplitude)
    {
        v.Unitize();
        return (v * Amplitude);
    }

    public static Vector3d RandomVector(int seed)
    {
        Vector3d r = Vector3d.Zero;
        Random rnd = new Random(seed);
        r = new Vector3d((rnd.NextDouble() * 2) - 1, (rnd.NextDouble() * 2) - 1, (rnd.NextDouble() * 2) - 1);
        return r;
    }

    public static Vector3d RandomVectorOnMesh(Mesh M, Point3d p, int seed)
    {
        Vector3d r = Vector3d.Zero;
        Vector3d MeshNormal = M.NormalAt(M.ClosestMeshPoint(p, 50));
        Plane VPlane = new Plane(p, MeshNormal);
        Random rnd = new Random(seed);
        double ang = rnd.NextDouble() * Math.PI * 2;
        r = VPlane.XAxis;
        r.Rotate(ang, VPlane.ZAxis);
        return r;
    }

    #endregion

    // </Custom additional code> 

    private List<string> __err = new List<string>(); //Do not modify this list directly.
    private List<string> __out = new List<string>(); //Do not modify this list directly.
    private RhinoDoc doc = RhinoDoc.ActiveDoc;       //Legacy field.
    private IGH_ActiveObject owner;                  //Legacy field.
    private int runCount;                            //Legacy field.

    public override void InvokeRunScript(IGH_Component owner, object rhinoDocument, int iteration, List<object> inputs, IGH_DataAccess DA)
    {
        //Prepare for a new run...
        //1. Reset lists
        this.__out.Clear();
        this.__err.Clear();

        this.Component = owner;
        this.Iteration = iteration;
        this.GrasshopperDocument = owner.OnPingDocument();
        this.RhinoDocument = rhinoDocument as Rhino.RhinoDoc;

        this.owner = this.Component;
        this.runCount = this.Iteration;
        this.doc = this.RhinoDocument;

        //2. Assign input parameters
        Mesh M = default(Mesh);
        if (inputs[0] != null)
        {
            M = (Mesh)(inputs[0]);
        }

        List<Point3d> P = null;
        if (inputs[1] != null)
        {
            P = GH_DirtyCaster.CastToList<Point3d>(inputs[1]);
        }


        //3. Declare output parameters
        object A = null;


        //4. Invoke RunScript
        RunScript(M, P, ref A);

        try
        {
            //5. Assign output parameters to component...
            if (A != null)
            {
                if (GH_Format.TreatAsCollection(A))
                {
                    IEnumerable __enum_A = (IEnumerable)(A);
                    DA.SetDataList(1, __enum_A);
                }
                else
                {
                    if (A is Grasshopper.Kernel.Data.IGH_DataTree)
                    {
                        //merge tree
                        DA.SetDataTree(1, (Grasshopper.Kernel.Data.IGH_DataTree)(A));
                    }
                    else
                    {
                        //assign direct
                        DA.SetData(1, A);
                    }
                }
            }
            else
            {
                DA.SetData(1, null);
            }

        }
        catch (Exception ex)
        {
            this.__err.Add(string.Format("Script exception: {0}", ex.Message));
        }
        finally
        {
            //Add errors and messages... 
            if (owner.Params.Output.Count > 0)
            {
                if (owner.Params.Output[0] is Grasshopper.Kernel.Parameters.Param_String)
                {
                    List<string> __errors_plus_messages = new List<string>();
                    if (this.__err != null) { __errors_plus_messages.AddRange(this.__err); }
                    if (this.__out != null) { __errors_plus_messages.AddRange(this.__out); }
                    if (__errors_plus_messages.Count > 0)
                        DA.SetDataList(0, __errors_plus_messages);
                }
            }
        }
    }
}