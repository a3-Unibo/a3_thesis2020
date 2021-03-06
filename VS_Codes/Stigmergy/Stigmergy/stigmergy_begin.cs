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

        2 main classes:
        . Particle System (and a class Particle)
        . Environment (Mesh with scalar and optional vector field)

        RunScript function structure:

        . check inputs & eventually bypass execution
        . initialise Environment
        . initialise Particle System
        . Update live variables
        . Update Particle System
        . Update Environment
        . Extract Output Geometries and Data

         */

        // return on null input
        if (M == null || P == null) return;


        // initialize on first execution or mesh change
        if (AS == null || M.Vertices.Count != ME.scalarField.Length)
        {
            // initialize particle system
            ME = new MeshEnvironment(M, dR, eR);
            AS = new AgentSystem(P, ME);
        }

        // restore initial values on reset
        if (reset)
        {
            // initialize particle system
            ME.RestoreScalarField();
            AS = new AgentSystem(P, ME);

        }

        if (go)
        {
            // update runtime variables
            AS.seekRadius = sR;
            AS.seekIntensity = sI;
            AS.sensDist = sensDist; // old value 3.0
            AS.sensAng = sensAng;
            ME.diffusionRate = (float)dR;
            ME.evaporationRate = (float)eR;

            // update simulation
            switch (method)
            {
                case 0:
                    AS.UpdateRTree();
                    break;
                case 1:
                    AS.UpdateJones();
                    break;
            }

            // update environment (diffusion + evaporation)
            ME.Update();

            // update component
            Component.ExpireSolution(true);
        }

        // . . . . . . .  extract geometries

        // particles positions and velocites
        AS.GetPointsVectors(out pts, out vecs);
        points = pts;
        vectors = vecs;

        // colored mesh
        outMCol = ME.GetColoredMesh();

        // debug mode
        if (debug)
        {
            switch (method)
            {
                case 0:
                    neigh = AS.GetNeighPts();
                    neighCol = AS.GetNeighBrightness();
                    break;
                case 1:
                    sensors = AS.SensorsOut();
                    break;
            }
        }

        // ............................................................................
        // </Custom code> 
    }

    // <Custom additional code> 
    // ............................................................................Global Variables
    public AgentSystem AS;
    public MeshEnvironment ME;
    public GH_Point[] pts;
    public GH_Vector[] vecs;

    // ............................................................................Classes

    // .......................................................... Agent System ........

    public class AgentSystem
    {
        // . . . . . . . . . . . . . . . . . . . . . . fields
        public List<Agent> Agents;
        public double seekRadius;
        public double seekIntensity;
        public double sensDist;
        public double sensAng;
        public MeshEnvironment ME;

        // . . . . . . . . . . . . . . . . . . . . . . constructor
        public AgentSystem(List<Point3d> positions, MeshEnvironment ME)
        {
            this.ME = ME;
            Agents = new List<Agent>();
            for (int i = 0; i < positions.Count; i++)
                Agents.Add(new Agent(this, positions[i], RandomVectorOnMesh(ME.M, positions[i], i) * 1.5));

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
             to avoid simultaneous write to the same location (surest way to get an error,
             even using so called "thread-safe" collections - especially if they are made
             of structures and not classes):
             . each particle saves its own closest point index as an internal field
             . the resulting contribution is written in a separate non-parallel loop
               that updates the scalar field.
             */
            Parallel.ForEach(Agents, ag =>

            {
                // find mesh closest point index
                ag.currentCPindex = ME.MeshPoints.ClosestPoint(ag.position);
                // Seek brightest point direction
                ag.SeekAngVis();

            });

            // write to scalar field
            foreach (Agent ag in Agents)
                ME.scalarField[ag.currentCPindex] = (float)Math.Min(1.0, ME.scalarField[ag.currentCPindex] + ag.PheroStrength);

            // update particles
            if (Agents.Count > 50)
            {
                Parallel.ForEach(Agents, ag =>
                {
                    ag.Update();
                });
            }
            else
                foreach (Agent ag in Agents)
                    ag.Update();

        }

        /// <summary>
        /// Update with RTree method (samples points within a sphere)
        /// </summary>
        /// <param name="M"></param>
        public void UpdateRTree()
        {
            foreach (Agent ag in Agents)
            {
                // clear particle neighbours list
                ag.neighbours.Clear();

                // perform RTree neighbour search
                ME.MeshRTree.Search(new Sphere(ME.M.ClosestPoint(ag.position + ag.velocity * sensDist), seekRadius),
                    (sender, args) => { ag.neighbours.Add(args.Id); });

                // Seek brightest point direction
                ag.SeekColor();

                //// update scalar field
                int cpInd = ME.MeshPoints.ClosestPoint(ag.position);
                ME.scalarField[cpInd] = (float)Math.Min(1.0, ME.scalarField[cpInd] + ag.PheroStrength);

            }

            if (Agents.Count > 50)
            {
                Parallel.ForEach(Agents, ag =>
               {
                   ag.Update();
               });
            }
            else
                foreach (Agent ag in Agents)
                    ag.Update();

        }

        #region output_data

        public void GetPointsVectors(out GH_Point[] pts, out GH_Vector[] vecs)
        {
            pts = new GH_Point[Agents.Count];
            vecs = new GH_Vector[Agents.Count];

            for (int i = 0; i < Agents.Count; i++)
            {
                pts[i] = new GH_Point(Agents[i].position);
                vecs[i] = new GH_Vector(Agents[i].velocity);
            }
        }

        public DataTree<GH_Point> GetNeighPts()
        {
            DataTree<GH_Point> ptsOut = new DataTree<GH_Point>();

            for (int i = 0; i < Agents.Count; i++)
                ptsOut.EnsurePath(new GH_Path(i));

            // parallelize for more than 50 particles
            if (Agents.Count > 50)
            {

                var x = Partitioner.Create(0, Agents.Count);

                Parallel.ForEach(x, (range, loopstate) =>
                {
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        for (int j = 0; j < Agents[i].neighbours.Count; j++)
                            ptsOut.Add(new GH_Point(ME.MeshPoints[Agents[i].neighbours[j]].Location), new GH_Path(i));
                    }
                });
            }
            else
            {
                for (int i = 0; i < Agents.Count; i++)
                {
                    for (int j = 0; j < Agents[i].neighbours.Count; j++)
                        ptsOut.Add(new GH_Point(ME.MeshPoints[Agents[i].neighbours[j]].Location), new GH_Path(i));
                }
            }

            return ptsOut;
        }

        public DataTree<GH_Number> GetNeighBrightness()
        {

            DataTree<GH_Number> briOut = new DataTree<GH_Number>();

            for (int i = 0; i < Agents.Count; i++)
                briOut.EnsurePath(new GH_Path(i));

            for (int i = 0; i < Agents.Count; i++)
            {
                for (int j = 0; j < Agents[i].neighbours.Count; j++)
                    briOut.Add(new GH_Number(ME.scalarField[Agents[i].neighbours[j]]), new GH_Path(i));
            }

            return briOut;
        }

        public DataTree<GH_Vector> SensorsOut()
        {
            DataTree<GH_Vector> sensOut = new DataTree<GH_Vector>();

            for (int i = 0; i < Agents.Count; i++)
            {
                GH_Path path = new GH_Path(i);
                sensOut.Add(new GH_Vector(Agents[i].sensorL), path);
                sensOut.Add(new GH_Vector(Agents[i].sensorC), path);
                sensOut.Add(new GH_Vector(Agents[i].sensorR), path);
            }
            return sensOut;
        }

        #endregion
    }

    // .......................................................... Particle ...............
    public class Agent
    {
        public Point3d position;
        public Vector3d velocity;
        public Vector3d acceleration;
        public Vector3d sensorC, sensorL, sensorR;
        public double visAng;
        public double MaxSpeed;
        public float PheroStrength;
        public int currentCPindex;
        public List<int> neighbours;
        public AgentSystem agSys;

        public Agent(AgentSystem agSys, Point3d position, Vector3d velocity)
        {
            this.agSys = agSys;
            this.position = position;
            this.velocity = velocity;
            acceleration = Vector3d.Zero;
            MaxSpeed = 1.5f;
            PheroStrength = 0.1f;
            visAng = agSys.sensAng;
            neighbours = new List<int>();
        }

        public void Update()
        {
            Move();
        }

        public void SeekAngVis()
        {

        }

        public void SeekColor()
        {
            acceleration = velocity;

            if (neighbours.Count > 0)
            {
                acceleration = Vector3d.Zero;
                Vector3d desired = Vector3d.Zero;
                double bri;
                double MaxBri = -1.0;
                int MaxInd = -1;

                // find neighbours index with maximum brightness
                for (int i = 0; i < neighbours.Count; i++)
                {
                    bri = agSys.ME.scalarField[neighbours[i]];
                    if (bri > MaxBri)
                    {
                        MaxBri = bri;
                        MaxInd = i;
                    }
                }

                desired += agSys.ME.MeshPoints[neighbours[MaxInd]].Location - position;

                //desired.Unitize();
                //desired *= MaxSpeed;
                //acc = (desired - vel) * pSystem.seekIntensity;

                acceleration = desired * agSys.seekIntensity;
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

        void Move()
        {
            velocity = velocity * 0.5 + acceleration * 0.5;
            //vel = vel + acc;
            if (velocity.Length > MaxSpeed)
            {
                velocity.Unitize();
                velocity *= MaxSpeed;
            }
            if (velocity.Length < 1)
            {
                velocity.Unitize();
            }
            position += velocity;
        }

    }

    // .......................................................... Mesh Environment .......

    public class MeshEnvironment
    {
        public Mesh M;
        public float diffusionRate;
        public float evaporationRate;
        public float[] scalarField;
        readonly float[] initVal;
        public PointCloud MeshPoints;
        public RTree MeshRTree;
        readonly int[][] NeighbourMap;
        readonly float[] DiffusionWeights;


        // M must be a colored Mesh
        public MeshEnvironment(Mesh M, double diffusionRate, double evaporationRate)
        {
            // assign external variables
            this.M = M;
            this.diffusionRate = (float)diffusionRate;
            this.evaporationRate = (float)evaporationRate;
            // scalar field structures
            scalarField = PopulateScalarField(M);
            initVal = PopulateScalarField(M);
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

        void EvaporateField()
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

        int[][] BuildNeighboursMap(Mesh M)
        {
            int[][] NeighbourMap = new int[M.Vertices.Count][];

            Parallel.For(0, M.Vertices.Count, i =>
            {
                NeighbourMap[i] = M.Vertices.GetConnectedVertices(i);
            });

            return NeighbourMap;

        }

        float[] CalculateDiffusionWeights()
        {
            float[] Weights = new float[NeighbourMap.Length];

            Parallel.For(0, NeighbourMap.Length, i =>
            {
                Weights[i] = 1 / (float)NeighbourMap[i].Length;
            });

            return Weights;
        }

        void Diffusion()
        {

        }

        RTree PopulateMeshRTree(Mesh M)
        {
            RTree rt = new RTree();

            for (int i = 0; i < M.Vertices.Count; i++)
            {
                rt.Insert((Point3d)M.Vertices[i], i);
            }

            return rt;
        }

        PointCloud PopulatePointCloud(Mesh M)
        {
            PointCloud mP = new PointCloud();

            for (int i = 0; i < M.Vertices.Count; i++)
            {
                mP.Add(M.Vertices[i], M.Normals[i]);
            }

            return mP;
        }

        #region output_data_methods

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

        #endregion
    }


    // ............................................................................Utilities

    #region Utilities

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