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
using System.Linq;
// </Custom using> 


/// <summary>
/// This class will be instantiated on demand by the Script component.
/// </summary>
public class Script_Instance2 : GH_ScriptInstance // Script_Instance2 - name changed to avoid conflicts while writing code
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
    private void RunScript(bool reset, bool go, bool debug, int steps, List<System.Object> Pl, double sR, List<Polyline> body, ref object iter, ref object Bodies, ref object TipPoints)
    {
        // <Custom code> 

        /*
          agent bodies interaction - C#
          code by Alessio Erioli - (c) Co-de-iT 2019
          rev. 2 - 2020
        */

        // . . . . . . . . . . . . . . . . . . . . . . . . return on null data
        if (Pl == null || body == null) return;

        // variables for data extraction
        DataTree<Polyline> outBodies = new DataTree<Polyline>();


        // . . . . . . . . . . . . . . . . . . . . . . . . initialize system
        if (reset || aBS == null)
        {
            // cast input into Planes list
            List<Plane> plList = Pl.Select(p => (Plane)p).ToList();

            aBS = new AgentBodySystem(plList, body, sR);
            iterationsCount = 0;
        }


        if (go)
        {
            for (int i = 0; i < steps; i++)
            {
                // update System
                aBS.Update();

                // update iteration counter
                iterationsCount++;
            }
            Component.ExpireSolution(true);
        }


        // . . . . . . . . . . . . . . . . . . . . . . . . extract geometries and data

        // unfortunately, DataTree structures do not work with parallel writing
        // body output must be done in a sequential fashion
        for (int i = 0; i < aBS.Agents.Length; i++)
            outBodies.AddRange(aBS.Agents[i].ExtractBody(), new GH_Path(i));

        iter = iterationsCount;
        Bodies = outBodies;

        if (debug)
        {
            // ConcurrentBag is a thread-safe collection that admits parallel writing
            // but the order is scrambled (in this case order doesn't matter)
            ConcurrentBag<GH_Point> tPoints = new ConcurrentBag<GH_Point>();

            Parallel.For(0, aBS.Agents.Length, i =>
            {
                for (int j = 0; j < aBS.Agents[i].body.Tips.Count; j++)
                    tPoints.Add(new GH_Point(aBS.Agents[i].body.Tips[j].pos));
            });

            TipPoints = tPoints;
        }


        // </Custom code> 
    }

    // <Custom additional code> 

    // global variables
    public AgentBodySystem aBS;
    public int iterationsCount;

    // .................................................................... classes 
    // .............................................................................
    // .............................................................................

    // ...................................................... Simulation ..........
    // ............................................................................
    // ............................................................................
    public class AgentBodySystem
    {
        public AgentBody[] Agents;
        readonly RTree agentsRTree;
        public double searchRadius;
        public double globalDeviation;
        public double deviationThreshold;

        public AgentBodySystem(List<Plane> Pl, List<Polyline> body, double searchRadius)
        {
            Agents = new AgentBody[Pl.Count];
            agentsRTree = new RTree();
            this.searchRadius = searchRadius;
            deviationThreshold = 0.01; // for future implementation (stop when global deviation under this threshold)

            // . . . . . . . . . . . . . . . . . . . . . . . . build agents array & RTree
            // since agents are fixed it can be built in advance
            for (int i = 0; i < Pl.Count; i++)
            {
                Agents[i] = new AgentBody(Pl[i], body);
                agentsRTree.Insert(Agents[i].agentPlane.Origin, i);
            }

            // . . . . . . . . . . . . . . . . . . . . . . . . . . . build neighbours map
            FindNeighbours();
            foreach (AgentBody ab in Agents)
                ab.FindTipNeighbour();
        }

        public void Update()
        {
            globalDeviation = 0.0; // for future implementation

            Parallel.ForEach(Agents, ab =>
            {
                ab.Update();
            });

            // call tips UpdatePosition and UpdateDirection
            Parallel.ForEach(Agents, ab =>
            {
                foreach (Tip t in ab.body.Tips) t.Update();
            });
        }

        void FindNeighbours()
        {
            foreach (AgentBody ab in Agents)
                agentsRTree.Search(new Sphere(ab.body.O, searchRadius), (sender, args) =>
                {
                    if (Agents[args.Id] != ab) ab.Neighbours.Add(Agents[args.Id]);
                });
        }
    }

    // ...................................................... Agent Body ..........
    // ............................................................................
    // ............................................................................

    public class AgentBody
    {
        public Plane agentPlane;
        public Body body;
        public List<AgentBody> Neighbours;
        public AgentBodySystem ABS;

        public AgentBody(Plane agentPlane, List<Polyline> polylines)
        {
            // define agent body and orient it
            this.agentPlane = agentPlane;
            body = new Body(polylines);
            body.Orient(Plane.WorldXY, this.agentPlane);

            Neighbours = new List<AgentBody>(); // this list will be filled later by the FindNeighbours function in AgentBodySystem
        }

        public void Update()
        {
            TipsCohesion();
        }

        public void FindTipNeighbour()
        {
            /*
             each tip finds it own designed neighbour. A neighbour AgentPlane's Tip is a candidate neighbour if:
             . neighbour's tip is within the angle of vision 
             (angle between the current tip dir and the vector connecting the tip with the nighbour's tip
             is within the angle of vision)
             OR
             . neighbour's direction is within the angle of vision
             (angle between the tip direction and the reverse of the nighbour's tip direction is within the angle of vision)

            check which angle is smaller, calculate the distance from the tip to the neighbour's tip or prev respectively
            and select the candidate closer to the tip (find candidate with minimum distance)

             */

            double dSqPos, dSqPrev;
            double anglePos, angleDir;
            double minD;
            foreach (Tip tip in body.Tips)
            {
                minD = double.MaxValue; // 1.0 - max Search Radius (squared)
                foreach (AgentBody nearAg in Neighbours)
                {
                    foreach (Tip neighTip in nearAg.body.Tips)
                    {
                        // calculate angles
                        anglePos = Vector3d.VectorAngle(tip.dir, neighTip.pos - tip.pos);
                        angleDir = Vector3d.VectorAngle(tip.dir, neighTip.prev - tip.pos);

                        if (anglePos < angleDir && anglePos < tip.angleOfVis)
                        {
                            dSqPos = tip.pos.DistanceToSquared(neighTip.pos);
                            if (dSqPos < minD)
                            {
                                minD = dSqPos;
                                tip.neighbour = neighTip;
                            }
                        }
                        else if (angleDir < tip.angleOfVis)
                        {
                            dSqPrev = tip.pos.DistanceToSquared(neighTip.prev);
                            if (dSqPrev < minD)
                            {
                                minD = dSqPrev;
                                tip.neighbour = neighTip;
                            }
                        }
                    }
                }
            }


        }

        public void TipsCohesion()
        {
            foreach (Tip t in body.Tips) t.ComputeDesiredToFixed();
        }

        public List<Polyline> ExtractBody()
        {
            body.Rebuild();
            return body.Arms;
        }

    }

    // ............................................................ Body ..........
    // ............................................................................
    // ............................................................................
    public class Body
    {
        public List<Polyline> Arms;
        public List<Tip> Tips;
        public Point3d O;
        public Body(List<Polyline> polylines)
        {
            /*
             ATTENTION - data type vs reference type copy - or - shallow vs deep copies

             structures (such as Point3d, Vector3d) are passed BY DATA, so, in this case:
             Point3d newPoint = new Point3d(OldPoint);
             newPoint is a true (deep) copy of OldPoint

             classes (all Geometry types, such as Polylines, Curves, Surfaces, BReps) are
             passed BY REFERENCE, so, in this case:
             List<Polyline> Arms = new List<Polyline>(polylines);
             Arms is a SHALLOW copy, so it references always the same PolyLine list.
             This means that any change in the PolyLines in Arms is passed also to the original list,
             since references aren't true copies but just pointers to another object. Any change to a
             reference propagates to the original and all other references.
              
             in order to make a deep copy look for a Duplicate() or Copy() method;
             it may be the case sometimes that one has to write a custom method to ensure copy
            */

            // this makes a deep copy - CORRECT
            Arms = new List<Polyline>();
            foreach (Polyline p in polylines)
                Arms.Add(p.Duplicate());

            O = Arms[0][0];

            Tips = new List<Tip>();
            Point3d prev;
            int armIndex = 0;
            foreach (Polyline arm in Arms)
            {
                if (arm.Count < 3) prev = O; else prev = arm[arm.Count - 2];
                Tips.Add(new Tip(this, arm[arm.Count - 1], prev, armIndex));
            }
        }

        public void Orient(Plane oldPlane, Plane newPlane)
        {
            var x = Transform.PlaneToPlane(oldPlane, newPlane);

            O.Transform(x);

            for (int i = 0; i < Arms.Count; i++)
            {
                Arms[i].Transform(x);
                Tips[i].pos = Arms[i][Arms[i].Count - 1];
                Tips[i].prev = Arms[i][Arms[i].Count - 2];
                Tips[i].dir = Tips[i].pos - Tips[i].prev;
            }

        }

        public void Rebuild()
        {
            for (int i = 0; i < Arms.Count; i++)
                Arms[i][Arms[i].Count - 1] = Tips[i].pos;
        }
    }

    // ............................................................. Tip ..........
    // ............................................................................
    // ............................................................................
    public class Tip
    {
        public Body body;
        public Point3d pos;
        public Point3d prev;
        public Vector3d dir;
        public Vector3d desired;
        public Tip neighbour;
        public int armIndex;
        public double maxLen;
        public double cR;
        public double cI;
        public double angleOfVis;

        public Tip(Body body, Point3d pos, Point3d prev, int armIndex, double cR, double cI, double angleOfVis)
        {
            this.body = body;
            this.pos = pos;
            this.prev = prev;
            this.armIndex = armIndex;
            dir = pos - prev;
            maxLen = dir.Length * 2.0;
            this.cR = cR;
            this.cI = cI;
            this.angleOfVis = angleOfVis;
            desired = Vector3d.Zero;
        }

        // constructor chaining example                                                                  cR     cI     angleOfVis
        public Tip(Body body, Point3d pos, Point3d prev, int armIndex) : this(body, pos, prev, armIndex, 1.0f, 0.05f, Math.PI * 0.3) { }


        public void ComputeDesiredToFixed()
        {
            // reset desired vector
            desired = Vector3d.Zero;
            // if neighbour tips list is zero or null return
            if (neighbour == null) return;

            desired = neighbour.prev - pos;

        }

        public void Update()
        {
            UpdatePosition();
            UpdateDirection();
        }

        public void UpdatePosition()
        {
            pos += desired * cI;
            // limit tip position to a max reach
            Vector3d toTip = pos - prev;
            if (toTip.Length > maxLen)
            {
                toTip = Limit(toTip, maxLen);
                pos = prev + toTip;
            }
        }

        public void UpdateDirection()
        {
            dir = pos - prev;
        }

    }

    // .................................................................. Utilities 
    // .............................................................................

    public static Vector3d Limit(Vector3d v, double length)
    {
        Vector3d rv = new Vector3d(v);
        double len = rv.Length;
        if (len > length)
        {
            rv.Unitize();
            rv *= length;
        }
        return rv;
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
        List<Plane> Pl = null;
        if (inputs[0] != null)
        {
            Pl = GH_DirtyCaster.CastToList<Plane>(inputs[0]);
        }
        List<int> id = null;
        if (inputs[1] != null)
        {
            id = GH_DirtyCaster.CastToList<int>(inputs[1]);
        }
        List<Polyline> body = null;
        if (inputs[2] != null)
        {
            body = GH_DirtyCaster.CastToList<Polyline>(inputs[2]);
        }


        //3. Declare output parameters
        object A = null;


        //4. Invoke RunScript
        RunScript(Pl, id, body, ref A);

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
