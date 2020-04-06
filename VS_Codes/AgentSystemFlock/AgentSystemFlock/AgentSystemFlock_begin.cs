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
using System.Threading;
using System.Threading.Tasks;
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
    private void RunScript(bool reset, bool go, int mode, List<Point3d> P, List<Vector3d> V, double nR, double coS, double alS, double seS, double seR, ref object Ap, ref object Av)
    {
        // <Custom code>
        GH_Point[] ptsOut;
        GH_Vector[] vecOut;

        if (reset || AgSys == null)
            AgSys = new AgentSystem(P, V);

        if (go)
        {
            AgSys.NeighborhoodRadius = nR;
            AgSys.CohesionStrength = coS;
            AgSys.AlignmentStrength = alS;
            AgSys.SeparationStrength = seS;
            AgSys.SeparationRadius = seR;

            if (AgSys.Agents.Count < 600)
                AgSys.Update();
            else
                switch (mode)
                {
                    case 0:
                        AgSys.UpdateParallel();
                        break;

                    case 1:
                        AgSys.UpdateRTree();
                        break;
                }

            Component.ExpireSolution(true);
        }

        AgSys.GetPtsVecs(out ptsOut, out vecOut);
        Ap = ptsOut;
        Av = vecOut;

        // </Custom code>
    }

    // <Custom additional code> 
    public AgentSystem AgSys;

    public class AgentSystem
    {
        // fields
        public List<Agent> Agents;
        public double NeighborhoodRadius;
        public double CohesionStrength;
        public double AlignmentStrength;
        public double SeparationStrength;
        public double SeparationRadius;
        public double MaxSpeed;
        public double BoundingBoxSize;
        public double ContainmentStrength;
        RTree PointsRTree;

        // constructor
        public AgentSystem(List<Point3d> P, List<Vector3d> V)
        {
            Agents = new List<Agent>();

            for (int i = 0; i < P.Count; i++) // i+=1
            {
                Agent ag = new Agent(P[i], V[i]);
                ag.Flock = this;

                Agents.Add(ag);
            }

            MaxSpeed = 0.3;
            BoundingBoxSize = 30.0;
            ContainmentStrength = 1.0;
        }

        // methods

        public void Update()
        {
            // compute agent velocity variation


            // update velocity and position

        }

        public void UpdateParallel()
        {
            /*
             <access> <T> <name> (<args>) {}

            (<args>) => {}

            (IEnumerable<T> name, <T> var =>{})
            (List<Agent> Agents, Agent ag =>{})
            (Agents, ag =>{})
             */
            //

            // compute agent velocity variation


            // update velocity and position
        }

        public void UpdateRTree()
        {
            // build RTree


            // find neighbours for each agent with RTree and compute agent velocity variation


            // update agent velocity and position

        }

        public List<Agent> FindNeighbours(Agent ag)
        {
            List<Agent> neighbours = new List<Agent>();

            foreach (Agent neighbour in Agents)
                if (neighbour != ag && neighbour.position.DistanceTo(ag.position) < NeighborhoodRadius)
                    neighbours.Add(neighbour);

            return neighbours;
        }

        public void ComputeAgent(Agent ag)
        {
            // find neighbours

            // compute velocity variation

        }

        public void GetPtsVecs(out GH_Point[] pts, out GH_Vector[] vecs)
        {
            pts = new GH_Point[Agents.Count];
            vecs = new GH_Vector[Agents.Count];

            for (int i = 0; i < Agents.Count; i++)
            {
                pts[i] = new GH_Point(Agents[i].position);
                vecs[i] = new GH_Vector(Agents[i].velocity);
            }
        }
    }


    public class Agent
    {
        // fields
        public Point3d position;
        public Vector3d velocity;
        public Vector3d desiredVelocity;
        public AgentSystem Flock;

        // constructor
        public Agent(Point3d position, Vector3d velocity)
        {
            this.position = position;
            this.velocity = velocity;
            desiredVelocity = this.velocity;
        }

        // methods

        public void ComputeDesiredVelocity(List<Agent> neighbours)
        {
            // initialize desired velocity

            // ------------------------------ CONTAINMENT -------------------------------
            Containment();
            // -------------------------------- FLOCKING --------------------------------
            if (neighbours.Count > 0)
            {
                // . . . . . . . . . . . . . COHESION BEHAVIOUR  . . . . . . . . . . . . .

                // find neighbours average


                // go there


                // update desired velocity
 

                // . . . . . . . . . . . . . ALIGNMENT BEHAVIOUR  . . . . . . . . . . . . .


                // . . . . . . . . . . . . . SEPARATION BEHAVIOUR  . . . . . . . . . . . . .

            }

        }

        public void Containment()
        {


        }


        public void UpdateVelocityAndPosition()
        {
            // steering


            // limit velocity to MaxSpeed


            // update position

        }
    }
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
        bool reset = default(bool);
        if (inputs[0] != null)
        {
            reset = (bool)(inputs[0]);
        }

        bool go = default(bool);
        if (inputs[1] != null)
        {
            go = (bool)(inputs[1]);
        }

        List<Point3d> P = null;
        if (inputs[2] != null)
        {
            P = GH_DirtyCaster.CastToList<Point3d>(inputs[2]);
        }
        List<Vector3d> V = null;
        if (inputs[3] != null)
        {
            V = GH_DirtyCaster.CastToList<Vector3d>(inputs[3]);
        }
        double nR = default(double);
        if (inputs[4] != null)
        {
            nR = (double)(inputs[4]);
        }

        double coS = default(double);
        if (inputs[5] != null)
        {
            coS = (double)(inputs[5]);
        }

        double alS = default(double);
        if (inputs[6] != null)
        {
            alS = (double)(inputs[6]);
        }

        double seS = default(double);
        if (inputs[7] != null)
        {
            seS = (double)(inputs[7]);
        }

        double seR = default(double);
        if (inputs[8] != null)
        {
            seR = (double)(inputs[8]);
        }



        //3. Declare output parameters
        object Ap = null;
        object Av = null;


        //4. Invoke RunScript
        RunScript(reset, go, P, V, nR, coS, alS, seS, seR, ref Ap, ref Av);

        try
        {
            //5. Assign output parameters to component...
            if (Ap != null)
            {
                if (GH_Format.TreatAsCollection(Ap))
                {
                    IEnumerable __enum_Ap = (IEnumerable)(Ap);
                    DA.SetDataList(1, __enum_Ap);
                }
                else
                {
                    if (Ap is Grasshopper.Kernel.Data.IGH_DataTree)
                    {
                        //merge tree
                        DA.SetDataTree(1, (Grasshopper.Kernel.Data.IGH_DataTree)(Ap));
                    }
                    else
                    {
                        //assign direct
                        DA.SetData(1, Ap);
                    }
                }
            }
            else
            {
                DA.SetData(1, null);
            }
            if (Av != null)
            {
                if (GH_Format.TreatAsCollection(Av))
                {
                    IEnumerable __enum_Av = (IEnumerable)(Av);
                    DA.SetDataList(2, __enum_Av);
                }
                else
                {
                    if (Av is Grasshopper.Kernel.Data.IGH_DataTree)
                    {
                        //merge tree
                        DA.SetDataTree(2, (Grasshopper.Kernel.Data.IGH_DataTree)(Av));
                    }
                    else
                    {
                        //assign direct
                        DA.SetData(2, Av);
                    }
                }
            }
            else
            {
                DA.SetData(2, null);
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