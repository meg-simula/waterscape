
import sys

from dolfin_utils.meshconvert import xml_writer 
import re 


def run(input, output): 
  """ 
  convert between amira and fenics implemented as state machine 
  0 = look for starters 
  1 = vertex 
  2 = cells 
  """
  state = 0 
  ifile = open(ifilename)
  ofile = open(ofilename, "w")

  vertices =[] 
  cells = []

  while 1: 
    line = ifile.readline()
    if line.startswith("BEGIN NODES"): 
      state = 1  
    elif line.startswith("END NODES"): 
      state = 0  
    elif line.startswith("BEGIN COMPONENTS"): 
      state = 2  
    elif line.startswith("END COMPONENTS"): 
      break
    elif state == 1:   
    # example: 
    #*node(1,-37.168381,50.233715,-12.345407,0,0,0,0,0)
      m = re.match(r"^\*node\((.*)\)$", line)
      if m: 
        index, x, y, z, n1, n2, n3, n4, n5  = m.group(1).split(",") 
        vertices.append((index,x,y,z))

    elif state == 2:   
    # example: 
#    *tetra4(1,1,38366,38367,9253,43441)
      m = re.match(r"^\*tetra4\((.*)\)$", line)
      if m: 
        index, n1, v1, v2, v3, v4  = m.group(1).split(",") 
        cells.append((index, v1, v2, v3, v4))
     
  xml_writer.write_header_mesh(ofile, "tetrahedron", 3) 

  xml_writer.write_header_vertices(ofile, len(vertices))
  counter = 0 
  for v in vertices: 
    if not counter == int(v[0])-1: print "Warning! Numbering assumption invalidated", counter, v[0]  
    xml_writer.write_vertex(ofile, int(v[0])-1, float(v[1]), float(v[2]), float(v[3])) 
    counter += 1 
  xml_writer.write_footer_vertices(ofile)

  xml_writer.write_header_cells(ofile, len(cells))
  counter=0
  for c in cells: 
    if not counter == int(c[0])-1: print "Warning! Numbering assumption invalidated", counter, c[0]  
    xml_writer.write_cell_tetrahedron(ofile, int(c[0])-1, int(c[1])-1, int(c[2])-1, int(c[3])-1, int(c[4])-1) 
    counter += 1 
  xml_writer.write_footer_cells(ofile)

  xml_writer.write_footer_mesh(ofile)
 
#   not handled: 
#    *component(1,"?",1,1)
#    *tria3(1,1,9253,38366,38367)
#
#    *component(2,"Inside",1,2)
#    *tetra4(1,1,38366,38367,9253,43441)
#      
#    if line.startswith("

if __name__ == "__main__": 
  ifilename=sys.argv[1]
  ofilename=sys.argv[2]


  run(ifilename, ofilename)

