using System;
using System.Drawing;
using Grasshopper.Kernel;

namespace AgentBodiesOnMesh
{
    public class AgentBodiesOnMeshInfo : GH_AssemblyInfo
    {
        public override string Name
        {
            get
            {
                return "AgentBodiesOnMesh";
            }
        }
        public override Bitmap Icon
        {
            get
            {
                //Return a 24x24 pixel bitmap to represent this GHA library.
                return null;
            }
        }
        public override string Description
        {
            get
            {
                //Return a short string describing the purpose of this GHA library.
                return "";
            }
        }
        public override Guid Id
        {
            get
            {
                return new Guid("f3031a4d-2c56-4e5e-88c1-67eb2934b77b");
            }
        }

        public override string AuthorName
        {
            get
            {
                //Return a string identifying you or your company.
                return "";
            }
        }
        public override string AuthorContact
        {
            get
            {
                //Return a string representing your preferred contact details.
                return "";
            }
        }
    }
}
