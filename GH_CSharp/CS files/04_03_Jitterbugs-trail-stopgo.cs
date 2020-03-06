using System;
using System.Collections;
using System.Collections.Generic;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;



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
  private void RunScript(List<Point3d> P, bool reset, bool go, ref object A, ref object trails)
  {
    
    if (reset || bugs.Count == 0 || P.Count != bugs.Count) initBugs(P);

    if (go)
    {
      // updates our JitterBugs
      foreach(JitterBug a in bugs) a.update();

      // expiring solution equals to a timer forcing the update
      Component.ExpireSolution(true);
    }
    
    // extract data for output
    for(int i = 0; i < bugs.Count; i++)
    {
      pts[i] = bugs[i].pos;
      trs[i] = bugs[i].trail;
    }

    A = pts;
    trails = trs;
  }

  // <Custom additional code> 
  
  // instantiate object from the class
  //JitterBug a;
  List<JitterBug> bugs = new List<JitterBug>();
  List<Point3d> pts = new List<Point3d>();
  List<Polyline> trs = new List<Polyline>();

  // whenever we want to use in a class something declared outside we should declare it as static
  //static Random rnd = new Random(); // a random class


  // initialization function
  public void initBugs(List<Point3d> P)
  {
    bugs.Clear();
    pts.Clear();
    trs.Clear();
    foreach(Point3d p in P)
    {
      bugs.Add(new JitterBug(p));
      pts.Add(p);
      trs.Add(new Polyline());
    }

  }


  // JitterBug class

  public class JitterBug
  {

    // static field
    static int nBugs = 0;

    // instance fields of the class
    public Point3d pos;
    public Vector3d vel;
    public Polyline trail;
    int id;
    Random rnd;

    // constructor(s)
    public JitterBug()
    {
      id = nBugs;
      rnd = new Random(id);
      pos = new Point3d(rnd.NextDouble() - 0.5, rnd.NextDouble() - 0.5, 0);
      vel = new Vector3d(rnd.NextDouble() - 0.5, rnd.NextDouble() - 0.5, 0) * 0.25;
      trail = new Polyline();
      nBugs++;
    }

    public JitterBug(Point3d pos){
      id = nBugs;
      rnd = new Random(id);
      this.pos = pos;
      vel = new Vector3d(rnd.NextDouble() - 0.5, rnd.NextDouble() - 0.5, 0) * 0.25;
      trail = new Polyline();
      nBugs++;
    }

    // methods

    public void update()
    {
      calcVel();
      move();
    }

    void calcVel()
    {
      vel.Rotate((rnd.NextDouble() - 0.5) * Math.PI * 0.6, Vector3d.ZAxis);
    }

    // check what changes if we invoke calcVel2
    void calcVel2()
    {
      vel = new Vector3d(rnd.NextDouble() - 0.5, rnd.NextDouble() - 0.5, 0);
    }


    void move()
    {
      trail.Add(pos);
      pos += vel;
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
    this. doc = this.RhinoDocument;

    //2. Assign input parameters
        List<Point3d> P = null;
    if (inputs[0] != null)
    {
      P = GH_DirtyCaster.CastToList<Point3d>(inputs[0]);
    }
    bool reset = default(bool);
    if (inputs[1] != null)
    {
      reset = (bool)(inputs[1]);
    }

    bool go = default(bool);
    if (inputs[2] != null)
    {
      go = (bool)(inputs[2]);
    }



    //3. Declare output parameters
      object A = null;
  object trails = null;


    //4. Invoke RunScript
    RunScript(P, reset, go, ref A, ref trails);
      
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
      if (trails != null)
      {
        if (GH_Format.TreatAsCollection(trails))
        {
          IEnumerable __enum_trails = (IEnumerable)(trails);
          DA.SetDataList(2, __enum_trails);
        }
        else
        {
          if (trails is Grasshopper.Kernel.Data.IGH_DataTree)
          {
            //merge tree
            DA.SetDataTree(2, (Grasshopper.Kernel.Data.IGH_DataTree)(trails));
          }
          else
          {
            //assign direct
            DA.SetData(2, trails);
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