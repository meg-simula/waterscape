import scipy.io
from dolfin import *
import numpy as np

def create_submesh_markers(markers, submesh, name):

    # creates the submesh markers
    submesh_markers = MeshFunction("size_t", submesh, markers.dim())

    # submesh to mesh mappings
    cell_map = submesh.data().array("parent_cell_indices", 3)
    #vertex_map = submesh.data().mesh_function("parent_vertex_indices")

    # iterate over the cells
    for cell in cells(submesh):
        i = cell.index()
        j = int(cell_map[i])
        submesh_markers.array()[i] = markers[j]

        # extracts the facets and the corresponding one of the full mesh
        #facets_sub = cell.entities(2)
        #facets = Cell(mesh, cell_map.array()[cell.index()]).entities(2)

        #for facet_sub in facets_sub :
        #    # mapped vertices
        #    vert_mapped = pfun.array()[Facet(submesh, facet_sub).entities(0)]
        #    # find the corresponding facet on the full mesh
        #    for facet in facets :
        #        # vertices
        #        vert = Facet(mesh, facet).entities(0)
        #        # intersection
        #        common = set(vert).intersection(snp.savetxt('test.out', x, fmt='%1.4e')et(vert_mapped))
        #        if len(common) == 3 :
        #            markers_submesh.array()[facet_sub] = markers.array()[facet]
        #            break
    np.save(name + "_cellmap", cell_map)


    return submesh_markers

def split_meshes():

    mesh = Mesh("bi18.xml")
    markers = MeshFunction("size_t", mesh, "bi18_sub.xml")
    #plot(markers, title="Cell markers")

    # Create mesh of the white and gray matter
    RGM = 22 
    LGM = 23
    RWM = 9
    LWM = 10 
    tissue = CellFunction("size_t", mesh, 0)
    tissue.array()[markers.array()==RGM] = 1
    tissue.array()[markers.array()==LGM] = 1
    tissue.array()[markers.array()==RWM] = 1
    tissue.array()[markers.array()==LWM] = 1
    whitegray = SubMesh(mesh, tissue, 1)

    whitetissue = CellFunction("size_t", mesh, 0)
    whitetissue.array()[markers.array()==RWM] = 1
    whitetissue.array()[markers.array()==LWM] = 1    

    white = SubMesh(mesh, whitetissue, 1)
  
    graytissue = CellFunction("size_t", mesh, 0)
    graytissue.array()[markers.array()==LGM] = 1
    graytissue.array()[markers.array()==RGM] = 1     

    gray = SubMesh(mesh, graytissue, 1)
    
    #white = SubMesh(mesh, markers, WHITE)
    #gray = SubMesh(mesh, markers, GRAY)

    # Map supermesh markers to submesh
    whitegray_markers = create_submesh_markers(markers, whitegray, "whitegray")
    white_markers = create_submesh_markers(markers, white, "white")
    gray_markers = create_submesh_markers(markers, gray, "gray")

    # Store submesh
    file = File("bi18_whitegray_mesh.xml.gz")
    file << whitegray
    file = File("bi18_white_mesh.xml.gz")
    file << white
    file = File("bi18_gray_mesh.xml.gz")
    file << gray
    #file = File("adult_mni_152_model_whitegray_mesh.xdmf")
    #file << whitegray

    # Store submarkers
    file = File("bi18_whitegray_markers.xml.gz")
    file << whitegray_markers
    file = File("bi18_white_markers.xml.gz")
    file << white_markers
    file = File("bi18_gray_markers.xml.gz")
    file << gray_markers
    #file = File("adult_mni_152_model_whitegray_markers.xdmf")
    #ile << whitegray_markers

    #plot(white, title="White matter")
    #plot(gray, title="Gray matter")
    #plot(whitegray, title="Mesh")
    #plot(whitegray_markers, title="Gray and white matter")
    #interactive()
def write_fibers():
    mesh = Mesh("bi18_whitegray_mesh.xml.gz")
    DG = VectorFunctionSpace(mesh, "DG", 0)
    fibers = Function(DG)
    cellmap = np.load("whitegray_cellmap.npy")
    lwt = np.loadtxt("bi18_labels7_all_fit_lwm_tetra.info")
    rwt = np.loadtxt("bi18_labels7_all_fit_rwm_tetra.info")
    lgt = np.loadtxt("bi18_labels7_all_fit_lgm_tetra.info")
    rgt = np.loadtxt("bi18_labels7_all_fit_rgm_tetra.info")

    rwe = open("bi18_rwm_eigvect_230000.asc", "r")
    lwe = open("bi18_lwm_eigvect_230000.asc", "r")
    lge = open("bi18_rgm_eigvect_230000.asc", "r")
    rge = open("bi18_rwm_eigvect_230000.asc", "r") 
    
    lw_eigvect = [[line.strip().split(" ")[:3] for line in lwe.readlines()[:]]]
    print lw_eigvect

    #for i in range(len(lwt)):
    #    fibers.vector()[in(lwt[i]-1)][0] = 
        
    #for cell in cells(mesh)
    
def convert_mat_to_xml(prefix, gdim=3, celltype="tetrahedron"):

    # Load matlab data from Colin27 or Adult MNI mesh
    if True:
        mat_contents = scipy.io.loadmat('MMC_Collins_Atlas_Mesh_Version_2L.mat')
        # Extract into separate arrays
        node = mat_contents["node"]
        elem = mat_contents["elem"]
        face = mat_contents["face"]

    else:
        mat_contents = scipy.io.loadmat('HeadVolumeMesh.mat')
        node =  mat_contents["HeadVolumeMesh"]["node"][0][0]
        elem =  mat_contents["HeadVolumeMesh"]["elem"][0][0]
        face =  mat_contents["HeadVolumeMesh"]["face"][0][0]

    num_vertices = node.size
    print "num_vertices = ", num_vertices
    num_cells = elem.size
    print "num_cells = ", num_cells

    # -- Write mesh

    # NB note the -1 offset here (matlab starts counting at 1, python
    # starts counting at 0)
    cells = "\n".join(["""      <%s index="%d" v0="%d" v1="%d" v2="%d" v3="%d"/>""" % (celltype, i, c[0]-1, c[1]-1, c[2]-1, c[3]-1)
                          for (i, c) in enumerate(elem)])
    vertices = "\n".join(["""      <vertex index="%d" x="%f" y="%f" z="%f"/>""" % (i, v[0], v[1], v[2])
                          for (i, v) in enumerate(node)])

    contents = """
<?xml version="1.0" encoding="UTF-8"?>

<dolfin xmlns:dolfin="http://fenicsproject.org">
  <mesh celltype="%s" dim="%d">
    <vertices size="%d">
    %s
    </vertices>

    <cells size="%d">
    %s
    </cells>
  </mesh>
</dolfin>
    """ % (celltype, gdim, num_vertices, vertices, num_cells, cells)

    # -- Domain markers
    N = gdim+1 # Which row contains the markers?
    markers = "\n".join(["""      <value cell_index="%d" local_entity="0" value="%d" />""" % (i, c[N])
                         for (i, c) in enumerate(elem)])

    markers_code = """
<?xml version="1.0"?>

<dolfin xmlns:dolfin="http://fenicsproject.org">
  <mesh_function>
    <mesh_value_collection name="m" type="uint" dim="%d" size="%d">
    %s
    </mesh_value_collection>
  </mesh_function>
</dolfin>
    """ % (gdim, num_cells, markers)

    filename1 = prefix + "_mesh.xml"
    print "Writing data to %s" % filename1
    file = open(filename1, 'w')
    file.write(contents)
    file.close()

    filename2 = prefix + "_markers.xml"
    print "Writing data to %s" % filename2
    file = open(filename2, 'w')
    file.write(markers_code)
    file.close()

    print "Success!"

    return filename1, filename2

if __name__ == "__main__":

    if False:
        prefix = "colin27"
        f1, f2 = convert_mat_to_xml(prefix)

        # Read mesh back in
        f1 = prefix + "_mesh.xml"
        mesh = Mesh(f1)

        # Read markers back in
        f2 = prefix + "_markers.xml"
        markers = MeshFunction("size_t", mesh, f2)
 
        # Store as h5 as well
        file = HDF5File(mpi_comm_world(), prefix + ".xdmf", 'w')
        file.write(mesh, "/mesh")
        file.write(markers, "/markers")
        file.close()

        # Look at mesh and markers
        plot(mesh, title=prefix)
        plot(markers, title=prefix)
        interactive()

    split_meshes()
    write_fibers()