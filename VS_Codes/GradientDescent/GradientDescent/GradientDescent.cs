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
using Rhino.Collections;
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
    private void RunScript(List<Point3d> P, Mesh M, double step, double eInd, double eMax, bool go, bool reset, ref object Me, ref object Pos, ref object Trails)
    {
        // <Custom code>

        // reset or initialization
        if (reset || erosion == null || P.Count != erosion.parts.Count) erosion = new ErosionSim(M, P);

        if (go)
        {
            // update live variables
            erosion.step = step;
            erosion.erosionIndex = eInd;
            erosion.erosionMax = eMax;

            // update System
            erosion.Update();

            // expire solution
            Component.ExpireSolution(true);
        }
        // data extraction for output
        pts = new Point3d[P.Count];
        trs = new Polyline[P.Count];

        for (int i = 0; i < erosion.parts.Count; i++)
        {
            pts[i] = erosion.parts[i].pos;
            if (erosion.parts[i].trail.IsValid) trs[i] = erosion.parts[i].trail;
            else trs[i] = null;
        }

        Me = erosion.extractMesh();
        Pos = pts;
        Trails = trs;

        // </Custom code>
    }

    // <Custom additional code> 

    // global variables
    public ErosionSim erosion;
    public Point3d[] pts;
    public Polyline[] trs;

    // simulation classes
    public class ErosionSim
    {
        // fields
        public List<Particle> parts;
        public Mesh M;
        public Mesh Morig;
        public Point3dList vertexList;
        public double step;
        public double erosionIndex;
        public double erosionMax;

        // constructor
        public ErosionSim(Mesh M, List<Point3d> P)
        {
            this.M = M;
            Morig = new Mesh();
            Morig.CopyFrom(M);
            vertexList = new Point3dList(M.Vertices.ToPoint3dArray());
            parts = new List<Particle>();
            foreach (Point3d p in P) parts.Add(new Particle(this, p));
        }

        // methods
        public void Update()
        {
            foreach (Particle pa in parts)
                if (pa.alive) pa.Update();

            for (int i = 0; i < vertexList.Count; i++)
                M.Vertices.SetVertex(i, vertexList[i]);

            M.RebuildNormals();
        }

        public Mesh extractMesh()
        {
            return M;
        }
    }

    public class Particle
    {
        // fields
        public Point3d pos;
        public Vector3d vel;
        public Polyline trail;
        public bool alive;
        public ErosionSim er;

        // constructor
        public Particle(ErosionSim er, Point3d pos)
        {
            this.er = er;
            this.pos = pos;
            vel = Vector3d.Zero;
            trail = new Polyline();
            trail.Add(new Point3d(pos));
            alive = true;
        }

        // methods

        public void Update()
        {
            CalcVel();
            Erode2();
            Move();
        }

        void Erode()
        {
            // find closest Mesh vertex to particle
            int vi = er.vertexList.ClosestIndex(pos);
            // if vertex displacement < erosionMax move vertex
            if (er.M.Vertices[vi].DistanceTo(er.Morig.Vertices[vi]) < er.erosionMax)
            {
                Vector3d disp = -1 * (Vector3d)er.M.Normals[vi] * er.erosionIndex;
                er.vertexList[vi] += disp;
            }
        }

        void Erode2()
        {
            // find closest Mesh vertex to particle
            int vi = er.vertexList.ClosestIndex(pos);
            // find vertex neighbours indexes
            int[] vNeigh = er.M.Vertices.GetConnectedVertices(vi);

            // if vertex displacement < erosionMax move vertex
            if (er.M.Vertices[vi].DistanceTo(er.Morig.Vertices[vi]) < er.erosionMax)
            {
                Vector3d disp = -1 * (Vector3d)er.M.Normals[vi] * er.erosionIndex;
                er.vertexList[vi] += disp;
                // apply displacement to nieghbours
                for (int i = 0; i < vNeigh.Length; i++)
                    er.vertexList[vNeigh[i]] += disp * 0.25;
            }
        }

        void CalcVel()
        {
            // get closest point to the Mesh
            MeshPoint mP = er.M.ClosestMeshPoint(pos, 0.0);
            // evaluate Mesh normal at that point
            Vector3d mN = er.M.NormalAt(mP);
            // find tangent to Mesh (cross product of normal and Z)
            vel = Vector3d.CrossProduct(mN, Vector3d.ZAxis);
            // rotate 90 degrees
            vel.Rotate(Math.PI * 0.5, mN); //
            // unitize and set to step value
            vel.Unitize();
            vel *= er.step;
        }

        void Move()
        {
            // compute new position
            Point3d newPos = pos + vel;
            newPos = er.M.ClosestPoint(newPos);

            // check if alive condition stands
            if (Math.Abs(newPos.Z - pos.Z) < 0.01) alive = false;
            else
            {
                pos = newPos;
                trail.Add(pos);
            }
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
        List<Point3d> P = null;
        if (inputs[0] != null)
        {
            P = GH_DirtyCaster.CastToList<Point3d>(inputs[0]);
        }
        Mesh M = default(Mesh);
        if (inputs[1] != null)
        {
            M = (Mesh)(inputs[1]);
        }

        double step = default(double);
        if (inputs[2] != null)
        {
            step = (double)(inputs[2]);
        }

        bool go = default(bool);
        if (inputs[3] != null)
        {
            go = (bool)(inputs[3]);
        }

        bool reset = default(bool);
        if (inputs[4] != null)
        {
            reset = (bool)(inputs[4]);
        }



        //3. Declare output parameters
        object Pos = null;
        object Trails = null;


        //4. Invoke RunScript
        RunScript(P, M, step, go, reset, ref Pos, ref Trails);

        try
        {
            //5. Assign output parameters to component...
            if (Pos != null)
            {
                if (GH_Format.TreatAsCollection(Pos))
                {
                    IEnumerable __enum_Pos = (IEnumerable)(Pos);
                    DA.SetDataList(1, __enum_Pos);
                }
                else
                {
                    if (Pos is Grasshopper.Kernel.Data.IGH_DataTree)
                    {
                        //merge tree
                        DA.SetDataTree(1, (Grasshopper.Kernel.Data.IGH_DataTree)(Pos));
                    }
                    else
                    {
                        //assign direct
                        DA.SetData(1, Pos);
                    }
                }
            }
            else
            {
                DA.SetData(1, null);
            }
            if (Trails != null)
            {
                if (GH_Format.TreatAsCollection(Trails))
                {
                    IEnumerable __enum_Trails = (IEnumerable)(Trails);
                    DA.SetDataList(2, __enum_Trails);
                }
                else
                {
                    if (Trails is Grasshopper.Kernel.Data.IGH_DataTree)
                    {
                        //merge tree
                        DA.SetDataTree(2, (Grasshopper.Kernel.Data.IGH_DataTree)(Trails));
                    }
                    else
                    {
                        //assign direct
                        DA.SetData(2, Trails);
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