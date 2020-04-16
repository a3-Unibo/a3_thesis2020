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
using System.Linq;
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
    private void RunScript(bool reset, bool go, bool debug, List<System.Object> P, List<System.Object> V, double curlNoiseScale, double curlNoiseTime, double simNoiseScale, double simNoiseTime, double isoValue, double isoValueTol, double planeRadius, double searchRadius, double cInt, double aInt, double sInt, double fieldInt, double isoVInt, double maxForce, ref object Planes, ref object neigh)
    {
        // <Custom code> 

        // align planes with field/neighbours - C#
        // code by Alessio Erioli - (c) Co-de-iT 2019 

        if (P == null || V == null || P.Count != V.Count) return;

        if (reset || aPS == null)
        {
            // Casting P and V to Point3d and Vector3d Lists
            List<Point3d> pList = new List<Point3d>();
            List<Vector3d> vList = new List<Vector3d>();
            for (int i = 0; i < P.Count; i++)
            {
                pList.Add((Point3d)P[i]);
                vList.Add((Vector3d)V[i]);
            }

            // passing the essential parameters to the new simulation
            aPS = new AgentPlaneSimulation(pList, vList);

        }

        if (go)
        {

            // update parameters
            aPS.curlNoiseScale = curlNoiseScale;
            aPS.curlNoiseTime = curlNoiseTime;
            aPS.simplexNoiseScale = (float)simNoiseScale;
            aPS.simplexNoiseTime = (float)simNoiseTime;
            aPS.isoValue = (float)isoValue;
            aPS.isoValueTol = (float)isoValueTol;
            aPS.planeRadius = planeRadius;
            aPS.searchRadius = searchRadius;
            aPS.cohesionIntensity = cInt;
            aPS.alignmentIntensity = aInt;
            aPS.separationIntensity = sInt;
            aPS.fieldIntensity = fieldInt;
            aPS.isoValIntensity = isoVInt;
            aPS.maxForce = maxForce;

            // run simulation
            aPS.Update();
            Component.ExpireSolution(true);
        }

        // extract output geometries
        Planes = aPS.ExtractPlanes();
        if (debug) neigh = aPS.ExtractNeighbours();
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
        public float simplexNoiseScale;
        public float simplexNoiseTime;
        public float isoValue;
        public float isoValueTol;
        public double planeRadius;
        public double searchRadius;
        public double cohesionIntensity;
        public double alignmentIntensity;
        public double separationIntensity;
        public double fieldIntensity;
        public double isoValIntensity;
        public double maxForce;
        public double sepDistSquared;

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

            // populate agent planes array
            Parallel.For(0, P.Length, i =>
            {
                agentPlanes[i] = new AgentPlane(P[i], V[i], this);
            });

        }

        public void Update()
        {
            sepDistSquared = 4 * planeRadius * planeRadius; // the square of 2*plane radius

            // . . . . . . . . . . environmental and stigmergic interactions


            // calculate field influence vector for agents

            Parallel.For(0, agentPlanes.Length, i =>
             {
                 // this resets the acceleration (or desired direction)
                 // in case a different update sequence is implemented
                 // remember to reset desired before calculating the new iteration
                 agentPlanes[i].ResetDesireds();
                 // computes alignment with Curl Noise
                 agentPlanes[i].ComputeCurlNoiseAlign();
                 agentPlanes[i].ComputeDesiredIsoValue();

             });

            // . . . . . . . . . . peer-to-peer interaction
            FlockRTree();

            // . . . . . . . . . . update position and direction
            UpdateAgentsDirection();

        }

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

                rTree.Search(new Sphere(agent.O, searchRadius), (sender, args) =>
                {
                    if (agentPlanes[args.Id] != agent)
                        agent.neighbours.Add(agentPlanes[args.Id]);
                });
            }

            // compute desired vector for each agent
            Parallel.ForEach(agentPlanes, agent =>
             {
                 agent.ComputeDesiredNeighbours();
             });

        }

        public void UpdateAgentsDirection()
        {
            // another way of implementing a parallel loop: parallel foreach with partitioner
            var part = Partitioner.Create(0, agentPlanes.Length);
            Parallel.ForEach(part, (range, loopstate) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                    agentPlanes[i].Update();
            });
        }


        public GH_Plane[] ExtractPlanes()
        {
            GH_Plane[] gP = new GH_Plane[agentPlanes.Length];
            Parallel.For(0, agentPlanes.Length, i =>
            {
                gP[i] = new GH_Plane(agentPlanes[i].ExtractPlane());
            });
            return gP;
        }

        public DataTree<GH_Point> ExtractNeighbours()
        {
            DataTree<GH_Point> neighPts = new DataTree<GH_Point>();
            for (int i = 0; i < agentPlanes.Length; i++)
                for (int j = 0; j < agentPlanes[i].neighbours.Count; j++)
                    neighPts.Add(new GH_Point(agentPlanes[i].neighbours[j].O), new GH_Path(i));

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
        public Vector3d X, Y;
        public Vector3d desDirX, desDirY; // desired directions - influenced by environmental field and neighbour agents' directions
        public Vector3d desPos; // desired position - influenced by environmental movement constraints (es. surface) and neighbour agents
        public AgentPlaneSimulation agentSim;
        public List<AgentPlane> neighbours;
        public float currentFieldVal;
        // search icosahedron for isovalue detection - common to all AgentPlanes and constant so it is declared static and readonly
        static readonly Vector3d[] SearchIcosahedron = {
            new Vector3d( 0, 0, 1), new Vector3d(0.894427, 0, 0.447214), new Vector3d(0.276393, 0.850651, 0.447214),
            new Vector3d( -0.723607, 0.525731, 0.447214), new Vector3d( -0.723607, -0.525731, 0.447214),
            new Vector3d( 0.276393, -0.850651, 0.447214), new Vector3d( 0.723607, 0.525731, -0.447214),
            new Vector3d( -0.276393, 0.850651, -0.447214), new Vector3d( -0.894427, 0, -0.447214),
            new Vector3d( -0.276393, -0.850651, -0.447214), new Vector3d( 0.723607, -0.525731, -0.447214), new Vector3d( 0, 0, -1)};

        // constructor
        public AgentPlane(Point3d O, Vector3d X, AgentPlaneSimulation agentSim)
        {
            this.O = O;
            X.Unitize();
            this.X = X;
            this.agentSim = agentSim;
            Vector3d oV = new Vector3d(O);
            oV.Unitize();
            Y = Vector3d.CrossProduct(X, oV);
            neighbours = new List<AgentPlane>();
        }
        // methods

        public void Update()
        {
            // update position
            O += desPos * agentSim.maxForce;
            // update X direction
            X = X * (1 - agentSim.maxForce) + agentSim.maxForce * desDirX;
            X.Unitize();
            // update Y direction
            Y = Y * (1 - agentSim.maxForce) + agentSim.maxForce * desDirY;
            Y.Unitize();
        }

        public void ResetDesireds()
        {
            desDirX = Vector3d.Zero;
            desDirY = Vector3d.Zero;
            desPos = Vector3d.Zero;
        }

        public void ComputeDesiredNeighbours()
        {

            // neighbours interaction
            if (neighbours.Count != 0)
            {
                // define flocking vectors
                Vector3d alignX = Vector3d.Zero;
                Vector3d alignY = Vector3d.Zero;
                Vector3d cohesion = Vector3d.Zero;
                Vector3d separation = Vector3d.Zero;
                Vector3d pointer;

                int sepCount = 0;
                double distanceSquared;
                double invNCount = 1.0 / neighbours.Count;

                // neighbours calculations
                foreach (AgentPlane neighbour in neighbours)
                {
                    alignX += neighbour.X;
                    alignY += neighbour.Y;
                    pointer = neighbour.O - O;
                    cohesion += pointer;
                    distanceSquared = O.DistanceToSquared(neighbour.O);
                    if (distanceSquared < agentSim.sepDistSquared)
                    {
                        // separation vector is bigger when closer to another agent
                        separation += -pointer * (agentSim.sepDistSquared - distanceSquared);
                        sepCount++;
                    }
                }

                // ............................ alignment behavior
                alignX *= invNCount;
                alignY *= invNCount;

                // updates desired direction
                // multiplies alignment vector by alignment intensity factor
                desDirX += agentSim.alignmentIntensity * alignX;
                desDirY += agentSim.alignmentIntensity * alignY;

                // ............................ cohesion behavior
                cohesion *= invNCount;

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

        public void ComputeDesiredIsoValue()
        {
            float fieldVal, diff, minDiff; // prepare variables

            currentFieldVal = Noise.Generate((float)(O.X * agentSim.simplexNoiseScale + agentSim.simplexNoiseTime),
                (float)(O.Y * agentSim.simplexNoiseScale + agentSim.simplexNoiseTime),
                (float)(O.Z * agentSim.simplexNoiseScale + agentSim.simplexNoiseTime));// field value at current positions

            minDiff = Math.Abs(currentFieldVal - agentSim.isoValue); // update current distance from iso value

            Point3d target = O;// if no value is within the difference range, the target is the current position

            if (minDiff > agentSim.isoValueTol)
            {
                for (int i = 0; i < SearchIcosahedron.Length; i++)
                {
                    Point3f sample = (Point3f)(O + SearchIcosahedron[i] * agentSim.planeRadius);
                    fieldVal = Noise.Generate(sample.X * agentSim.simplexNoiseScale + agentSim.simplexNoiseTime,
                        sample.Y * agentSim.simplexNoiseScale + agentSim.simplexNoiseTime,
                        sample.Z * agentSim.simplexNoiseScale + agentSim.simplexNoiseTime);
                    diff = Math.Abs(fieldVal - agentSim.isoValue);
                    if (diff < minDiff)
                    {
                        minDiff = diff;
                        target = sample;
                    }
                }
            }

            desPos += (target - O) * agentSim.isoValIntensity;
        }

        public void ComputeCurlNoiseAlign()
        {
            // curl noise vector
            Vector3d curlNoiseVector;
            curlNoiseVector = CurlNoise.CurlNoiseVector(O, agentSim.curlNoiseScale, agentSim.curlNoiseTime, true);
            curlNoiseVector.Unitize();
            ComputeFieldAlign(curlNoiseVector);
        }

        /// <summary>
        /// This is a more generic function to allow the agent plane system to interact with more than one vector field.
        /// </summary>
        /// <param name="fieldVector">The Vector the agent interacts with</param>
        public void ComputeFieldAlign(Vector3d fieldVector)
        {
            desDirX += fieldVector * agentSim.fieldIntensity;
        }


        public Plane ExtractPlane()
        {
            //Vector3d oV = new Vector3d(O);
            //oV.Unitize();
            //Y = Vector3d.CrossProduct(X, oV);
            return new Plane(O, X, Y);
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