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

using System.Threading.Tasks;
using System.Collections.Concurrent;
using System.Drawing;
using Noises;
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
    private void RunScript(bool reset, bool go, bool debug, List<Point3d> P, List<Vector3d> V, double curlNoiseScale, double curlNoiseTime, double perlNoiseScale, double perlNoiseTime, double isoValue, double isoValueTol, double planeRadius, double searchRadius, double cInt, double aInt, double sInt, double fieldInt, double isoVInt, double maxForce, ref object Planes, ref object neigh)
    {
        // <Custom code> 

        // align planes with field/neighbours - C#
        // code by Alessio Erioli - (c) Co-de-iT 2019 


        if (P == null || P.Count == 0) return;

        if (reset || aPS == null)
        {
            // passing the essential parameters to the new simulation
            aPS = new AgentPlaneSimulation(P, V);

        }

        // update parameters
        aPS.curlNoiseScale = curlNoiseScale;
        aPS.curlNoiseTime = curlNoiseTime;
        aPS.perlNoiseScale = perlNoiseScale;
        aPS.perlNoiseTime = perlNoiseTime;
        aPS.isoValue = isoValue;
        aPS.isoValueTol = isoValueTol;
        aPS.planeRadius = planeRadius;
        aPS.searchRadius = searchRadius;
        aPS.cohesionIntensity = cInt;
        aPS.alignmentIntensity = aInt;
        aPS.separationIntensity = sInt;
        aPS.fieldIntensity = fieldInt;
        aPS.isoValIntensity = isoVInt;
        aPS.maxForce = maxForce;

        if (go)
        {
            // run simulation
            aPS.Run();
            Component.ExpireSolution(true);
        }

        // extract output geometries
        Planes = aPS.OutPlanes();
        if (debug) neigh = aPS.NeighOut();
        // </Custom code> 
    }

    // <Custom additional code> 

    // ............................................................ global variables 
    // global variables
    public AgentPlaneSimulation aPS;

    // .................................................................... classes 
    // .............................................................................
    // .............................................................................



    // ...................................................... Simulation ..........
    // ............................................................................
    // ............................................................................

    public class AgentPlaneSimulation
    {
        // ..........................    fields

        public double curlNoiseScale;
        public double curlNoiseTime;
        public double perlNoiseScale;
        public double perlNoiseTime;
        public double isoValue;
        public double isoValueTol;
        public double planeRadius;
        public double searchRadius;
        public double cohesionIntensity;
        public double alignmentIntensity;
        public double separationIntensity;
        public double fieldIntensity;
        public double isoValIntensity;
        public double maxForce;
        
        public AgentPlane[] agentPlanes;

        // ..........................    constructor

        public AgentPlaneSimulation(List<Point3d> P, List<Vector3d> V)
        {
            // build agent planes array
            Build(P.ToArray(), V.ToArray());
        }

        // ..........................    methods

        public void Build(Point3d[] P, Vector3d[] V)
        {
            // build agent planes array
            agentPlanes = new AgentPlane[P.Length];

            if (V == null || V.Length < P.Length)
                V = RandomVectors(P.Length);

            Parallel.For(0, P.Length, i =>
            {
                agentPlanes[i] = new AgentPlane(P[i], V[i], this);
            });

        }

        public void Run()
        {
            UpdateAgents();
        }

        public void UpdateAgents()
        {

            // . . . . . . . . . . environmental and stigmergic interactions


            // calculate field influence vector for agents

            Parallel.For(0, agentPlanes.Length, i =>
             {
                 // curl noise vector
                 Vector3d cNV;
                 cNV = CurlNoise.CurlNoiseVector(agentPlanes[i].O, curlNoiseScale, curlNoiseTime, true);
                 cNV.Unitize();

                 // this resets the acceleration (or desired direction)
                 // in case a different update sequence is implemented
                 // remember to reset desired before calculating the new iteration
                 agentPlanes[i].ResetDesireds();
                 agentPlanes[i].ComputeFieldAlign(cNV);
                 agentPlanes[i].ComputeDesiredIsoValue((float)isoValue, perlNoiseScale, perlNoiseTime);

             });

            // . . . . . . . . . . peer-to-peer interaction
            FlockRTree();

            // . . . . . . . . . . update position and direction
            UpdateAgentsDirection();

        }

        //public void UpdateAgents_OLD()
        //{
        // calculate curl noise vector for agents

        //Vector3d[] curlNoiseVector = new Vector3d[agentPlanes.Length];

        // parallel loop version 1
        //var p = Partitioner.Create(0, agentPlanes.Length);

        //Parallel.ForEach(p, (range, loopState) =>
        //{
        //    for (int i = range.Item1; i < range.Item2; i++)
        //    {
        //        CurlNoise(agentPlanes[i].O, cNS, cNT, true, out curlNoiseVector[i]);
        //        curlNoiseVector[i].Unitize();
        //        agentPlanes[i].ResetDesired();
        //        agentPlanes[i].AlignWithField(curlNoiseVector[i], fI);
        //    }
        //});


        // parallel loop version 2
        //Parallel.For(0, agentPlanes.Length, i =>
        //{
        //    CurlNoise(agentPlanes[i].O, cNS, cNT, true, out curlNoiseVector[i]);
        //    curlNoiseVector[i].Unitize();
        //    agentPlanes[i].Align(curlNoiseVector[i], mF);
        //});
        //}

        public void FlockRTree()
        {

            // declare RTree
            RTree rTree = new RTree();

            // populate RTree
            for (int i = 0; i < agentPlanes.Length; i++)
                rTree.Insert(agentPlanes[i].O, i);

            // find neighbours for each agent plane
            foreach (AgentPlane agent in agentPlanes)
            {
                agent.neighbours.Clear();

                EventHandler<RTreeEventArgs> rTreeCallback =
                    (object sender, RTreeEventArgs args) =>
                    {
                        if (agentPlanes[args.Id] != agent)
                            agent.neighbours.Add(agentPlanes[args.Id]);
                    };

                rTree.Search(new Sphere(agent.O, searchRadius), rTreeCallback);
            }

            // compute desired vector for each agent
            Parallel.ForEach(agentPlanes, agent =>
             {
                 agent.ComputeDesiredNeighbours();
             });

        }

        public void UpdateAgentsDirection()
        {
            var part = Partitioner.Create(0, agentPlanes.Length);
            Parallel.ForEach(part, (range, loopstate) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    // update at max Force
                    agentPlanes[i].Update(maxForce);
                }
            });
        }

        /// <summary>
        /// Output Planes
        /// </summary>
        /// <returns></returns>
        public GH_Plane[] OutPlanes()
        {
            GH_Plane[] gP = new GH_Plane[agentPlanes.Length];
            Parallel.For(0, agentPlanes.Length, i =>
            {
                gP[i] = new GH_Plane(agentPlanes[i].PlaneOut());

            });
            return gP;
        }

        public DataTree<GH_Point> NeighOut()
        {
            DataTree<GH_Point> neighPts = new DataTree<GH_Point>();
            for (int i = 0; i < agentPlanes.Length; i++)
                for (int j = 0; j < agentPlanes[i].neighbours.Count; j++)
                {
                    neighPts.Add(new GH_Point(agentPlanes[i].neighbours[j].O), new GH_Path(i));

                }
            return neighPts;
        }

    }


    // ........................................................... Agent ..........
    // ............................................................................
    // ............................................................................
    public class AgentPlane
    {

        // fields
        public Point3d O;
        public Vector3d X, Y, Z;
        public Vector3d desDir; // desired direction - influenced by environmental field and neighbour agents' directions
        public Vector3d desPos; // desired position - influenced by environmental movement constraints (es. surface) and neighbour agents
        public AgentPlaneSimulation agentSim;
        public List<AgentPlane> neighbours;
        public float currentFieldVal;
        public Point3d[] SearchIcosahedron =
        {
            new Point3d( 0, 0, 1),
            new Point3d(0.894427, 0, 0.447214),
            new Point3d(0.276393, 0.850651, 0.447214),
            new Point3d( -0.723607, 0.525731, 0.447214),
            new Point3d( -0.723607, -0.525731, 0.447214),
            new Point3d( 0.276393, -0.850651, 0.447214),
            new Point3d( 0.723607, 0.525731, -0.447214),
            new Point3d( -0.276393, 0.850651, -0.447214),
            new Point3d( -0.894427, 0, -0.447214),
            new Point3d( -0.276393, -0.850651, -0.447214),
            new Point3d( 0.723607, -0.525731, -0.447214),
            new Point3d( 0, 0, -1)};

        // constructor
        public AgentPlane(Point3d O, Vector3d dirX, AgentPlaneSimulation agentSim)
        {
            this.O = O;
            dirX.Unitize();
            this.X = dirX;
            this.agentSim = agentSim;
            neighbours = new List<AgentPlane>();
        }
        // methods

        public void Update(double maxForce)
        {
            // update position
            O += desPos * maxForce;
            // update direction
            X = X * (1 - maxForce) + maxForce * desDir;
            X.Unitize();
        }

        public void ResetDesireds()
        {
            desDir = Vector3d.Zero;
            desPos = Vector3d.Zero;
        }

        public void ComputeDesiredNeighbours()
        {

            // neighbours interaction
            if (neighbours.Count != 0)
            {
                // define flocking vectors
                Vector3d align = Vector3d.Zero;
                Vector3d cohesion = Vector3d.Zero;
                Vector3d separation = Vector3d.Zero;
                Vector3d pointer = Vector3d.Zero;
                int sepCount = 0;
                double sepDistSq = 4 * agentSim.planeRadius * agentSim.planeRadius; // the square of 2*plane radius
                double distSq;

                // neighbours calculations
                foreach (AgentPlane neighbour in neighbours)
                {
                    align += neighbour.X;
                    pointer = (Vector3d)Point3d.Subtract(neighbour.O, O);
                    cohesion += pointer;
                    distSq = O.DistanceToSquared(neighbour.O);
                    if (distSq < sepDistSq)
                    {
                        // separation vector is bigger when closer to another agent
                        separation += -pointer * (sepDistSq - distSq);
                        sepCount++;
                    }
                }

                // ............................ alignment behavior
                align /= neighbours.Count;

                // updates desired direction
                // multiplies alignment vector by alignment intensity factor
                desDir += agentSim.alignmentIntensity * align;

                // ............................ cohesion behavior
                cohesion /= neighbours.Count;

                // updates desired position
                // multiplies cohesion vector by cohesion intensity factor
                desPos += agentSim.cohesionIntensity * cohesion;


                // ............................ separation behavior
                if (sepCount > 0)
                    separation /= sepCount;

                // updates desired position
                // multiplies separation vector by separation intensity factor
                desPos += agentSim.separationIntensity * separation;

            }

        }

        public void ComputeDesiredIsoValue(float isoValue, double perlNoiseScale, double perlNoiseTime)
        {
            float fieldVal, diff, minDiff; // prepare variables

            currentFieldVal = Noise.Generate((float)(O.X * perlNoiseScale + perlNoiseTime),
                (float)(O.Y * perlNoiseScale + perlNoiseTime),
                (float)(O.Z * perlNoiseScale + perlNoiseTime));// field value at current positions

            minDiff = Math.Abs(currentFieldVal - isoValue); // update current distance from iso value

            Point3d target = O;// if no value is within the difference range, the target is the actual position

            if (minDiff > agentSim.isoValueTol)
            {
                for (int i = 0; i < SearchIcosahedron.Length; i++)
                {
                    Point3d sample = O + SearchIcosahedron[i] * agentSim.planeRadius;
                    fieldVal = Noise.Generate((float)(sample.X * perlNoiseScale + perlNoiseTime),
                        (float)(sample.Y * perlNoiseScale + perlNoiseTime),
                        (float)(sample.Z * perlNoiseScale + perlNoiseTime));
                    diff = Math.Abs(fieldVal - isoValue);
                    if (diff < minDiff)
                    {
                        minDiff = diff;
                        target = sample;
                    }
                }
            }

            desPos += (target - O) * agentSim.isoValIntensity;
        }

        public void ComputeFieldAlign(Vector3d fieldVector)
        {
            desDir += fieldVector * agentSim.fieldIntensity;
        }

        public void AlignWithVector(Vector3d vector, double intensity)
        {
            X = X * (1 - intensity) + vector * intensity;
            X.Unitize();
        }

        public Plane PlaneOut()
        {
            Vector3d oV = new Vector3d(O);
            oV.Unitize();
            this.Y = Vector3d.CrossProduct(X, oV);
            return new Plane(O, X, Y);
        }

        public GH_Mesh MeshOut(double rad)
        {
            Plane pl = PlaneOut();
            Transform x = Transform.PlaneToPlane(Plane.WorldXY, pl);

            Mesh m = new Mesh();
            Point3d[] verts = new Point3d[] { new Point3d(-1, -1, 0), new Point3d(1, -1, 0), new Point3d(1, 1, 0), new Point3d(-1, 1, 0) };
            m.Vertices.AddVertices(verts);
            m.Faces.AddFace(0, 1, 2, 3);
            m.Scale(rad);
            m.Transform(x);
            m.Normals.ComputeNormals();

            GH_Mesh m1 = new GH_Mesh(m);
            return m1;
        }
    }

    // .................................................................. Utilities 
    #region Utilities

    public float Remap(float val, float from1, float to1, float from2, float to2)
    {
        return (val - from1) / (to1 - from1) * (to2 - from2) + from2;
    }

    public GH_Colour Vec2Col(Vector3d v)
    {
        v.Unitize();
        int red = (int)Math.Floor(Remap((float)v.X, -1f, 1f, 0f, 255f));
        int green = (int)Math.Floor(Remap((float)v.Y, -1f, 1f, 0f, 255f));
        int blue = (int)Math.Floor(Remap((float)v.Z, -1f, 1f, 0f, 255f));
        return new GH_Colour(Color.FromArgb(red, green, blue));
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

    public static Vector3d[] RandomVectors(int length)
    {
        Vector3d[] randomVectors = new Vector3d[length];
        Parallel.For(0, length, i =>
        {
            randomVectors[i] = RandomVector(i);
        });

        return randomVectors;
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
        List<Point3d> P = null;
        if (inputs[0] != null)
        {
            P = GH_DirtyCaster.CastToList<Point3d>(inputs[0]);
        }
        List<Vector3d> V = null;
        if (inputs[1] != null)
        {
            V = GH_DirtyCaster.CastToList<Vector3d>(inputs[1]);
        }
        List<Vector3d> vF = null;
        if (inputs[2] != null)
        {
            vF = GH_DirtyCaster.CastToList<Vector3d>(inputs[2]);
        }
        double pR = default(double);
        if (inputs[3] != null)
        {
            pR = (double)(inputs[3]);
        }

        double cI = default(double);
        if (inputs[4] != null)
        {
            cI = (double)(inputs[4]);
        }

        double aI = default(double);
        if (inputs[5] != null)
        {
            aI = (double)(inputs[5]);
        }

        double sI = default(double);
        if (inputs[6] != null)
        {
            sI = (double)(inputs[6]);
        }

        object fI = default(object);
        if (inputs[7] != null)
        {
            fI = (object)(inputs[7]);
        }

        double mF = default(double);
        if (inputs[8] != null)
        {
            mF = (double)(inputs[8]);
        }



        //3. Declare output parameters
        object Planes = null;


        //4. Invoke RunScript
        RunScript(P, V, vF, pR, cI, aI, sI, fI, mF, ref Planes);

        try
        {
            //5. Assign output parameters to component...
            if (Planes != null)
            {
                if (GH_Format.TreatAsCollection(Planes))
                {
                    IEnumerable __enum_Planes = (IEnumerable)(Planes);
                    DA.SetDataList(1, __enum_Planes);
                }
                else
                {
                    if (Planes is Grasshopper.Kernel.Data.IGH_DataTree)
                    {
                        //merge tree
                        DA.SetDataTree(1, (Grasshopper.Kernel.Data.IGH_DataTree)(Planes));
                    }
                    else
                    {
                        //assign direct
                        DA.SetData(1, Planes);
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