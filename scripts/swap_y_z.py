import bpy

def swap_y_z_coordinates():
    # Get the active object
    obj = bpy.context.active_object

    # Check if an object is selected and is a mesh
    if obj and obj.type == 'MESH':
        # Get the mesh data
        mesh = obj.data
        
        # Iterate through vertices
        for vertex in mesh.vertices:
            # Store original coordinates
            old_y = vertex.co.y
            old_z = vertex.co.z
            
            # Swap x and y
            vertex.co.y = old_z
            vertex.co.z = old_y
            
        print(f"Swapped Y/Z for: {obj.name}")
    else:
        print("No mesh object selected.")

# Run the function
swap_y_z_coordinates()
