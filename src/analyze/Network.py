import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
import os

# matplotlib.use('TkAgg')


def main():

    analyze_dir = "analyze"
    outputFile = "output/analyze/components2.txt"

    out = open(outputFile, "w")
    for filename in os.listdir(analyze_dir):
        g = nx.Graph()

        # for i in range(0,2665):
        #     g.add_node(i)

        with open(os.path.join(analyze_dir, filename)) as f:
            for line in f:
                splitLine = line.split()
                edge1 = int(splitLine[0])
                edge2 = int(splitLine[1])
                weight = int(splitLine[2])
                if edge1 != edge2:
                    g.add_edge(edge1, edge2, weight=weight)
                # g.add_edge(edge1, edge2, weight=weight)

        out.write("Cell: {}\n".format(filename))
        out.write("# of edges: {}\n".format(g.number_of_edges()))
        out.write("# of nodes: {}\n".format(g.number_of_nodes()))
        out.write("Connected: {}\n".format(nx.is_connected(g)))
        out.write("Connected Components: {}\n".format(nx.number_connected_components(g)))
        out.write("\n")

        nx.draw(g)
        plt.show()

        # print("Cell: {}".format(filename))
        # print("# of edges: {}".format(g.number_of_edges()))
        # print("# of nodes: {}".format(g.number_of_nodes()))
        # print("Connected: {}".format(nx.is_connected(g)))
        # print("Connected Components: {}".format(nx.number_connected_components(g)))
        # print("\n")

    out.close()

if __name__ == '__main__':
    main()
