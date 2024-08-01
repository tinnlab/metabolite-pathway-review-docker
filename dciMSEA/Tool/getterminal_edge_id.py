

def getterminal_edge_id(terminal_edges_name,namelist):
    terminal_edges = []
    for i in terminal_edges_name:
        a = [j for j, x in enumerate(namelist) if x in i[0]]
        b = [j for j, x in enumerate(namelist) if x in i[1]]
        id = a + b
        if len(id) == 2:  # Ensure both indices are found
            terminal_edges.append(tuple(id))


    return  terminal_edges