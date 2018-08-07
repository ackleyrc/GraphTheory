using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class Graph
{
    HashSet<Node> nodes = new HashSet<Node>();
    HashSet<Edge> edges = new HashSet<Edge>();

    NodeEdgePairs edgesFromNode = new NodeEdgePairs();
    NodeEdgePairs edgesToNode = new NodeEdgePairs();

    public void AddUndirectedEdge(Edge newEdge)
    {
        AddDirectedEdge(newEdge);
        AddDirectedEdge(newEdge.GetReversed());
    }

    public void AddDirectedEdge(Edge newEdge)
    {
        if (edges.Add(newEdge))
        {
            edgesFromNode.Add(newEdge.GetFirstNode(), newEdge);
            edgesToNode.Add(newEdge.GetLastNode(), newEdge);
        }
    }

    /// <summary>
    /// Finds the shortest series of edges that continuously connect from the startNode to the endNode. 
    /// Optionally accepts a predicate function filtering out edges that do not meet the provided criteria. 
    /// If no path exists, the returned list will be empty. This implementation is based on the A* algorithm.
    /// </summary>
    /// <param name="startNode"></param>
    /// <param name="endNode"></param>
    /// <param name="criteria"></param>
    /// <returns></returns>
    public List<Node> GetShortestPath(Node startNode, Node endNode, NodeEdgePairs.EdgeFilterCriteria criteria = null)
    {                 
        List<Node> shortestPath = new List<Node>();

        HashSet<Node> openNodes = new HashSet<Node>(); // set of nodes to be evaluated
        HashSet<Node> closedNodes = new HashSet<Node>(); // set of nodes already evaluated

        Dictionary<Node, Node> cameFrom = new Dictionary<Node, Node>(); // associates a key node with the node it came from
        Dictionary<Node, float> gScore = new Dictionary<Node, float>(); // cost of getting to node from the startNode (inf if not included)
        Dictionary<Node, float> fScore = new Dictionary<Node, float>(); // total cost of getting from start to end via a given node (partly known, partly heuristic)

        openNodes.Add(startNode);
        gScore.Add(startNode, 0f); // cost from start to start is zero
        fScore.Add(startNode, (endNode.GetXYZ() - startNode.GetXYZ()).magnitude); // for start node, total cost is all heuristic (i.e unweighted distance to end)

        while (openNodes.Count > 0)
        {
            // Get the node in the open set that has the current lowest fScore
            Node currentNode = openNodes.First();
            foreach (Node node in openNodes)
            {
                if (fScore.ContainsKey(node) && fScore[node] < fScore[currentNode])
                {
                    currentNode = node;
                }
            }

            // if the current node is the target/goal endNode, then we've reached the end and can return the shortest path
            if (currentNode == endNode)
            {
                shortestPath.Add(currentNode);
                while (cameFrom.ContainsKey(currentNode))
                {
                    currentNode = cameFrom[currentNode];
                    shortestPath.Add(currentNode); // adds nodes from finish to start
                }
                return shortestPath.AsEnumerable().Reverse().ToList(); // reverse nodes to order from start to finish
            }

            openNodes.Remove(currentNode);
            closedNodes.Add(currentNode);

            foreach (Edge edgeFromNode in edgesFromNode.GetEdges(currentNode, criteria))
            {
                Node otherNode = edgeFromNode.GetOtherNode(currentNode);

                if (closedNodes.Contains(otherNode)) continue; // ignore if already evaluated
                
                float tentativeGScore = gScore[currentNode] + edgeFromNode.GetPathLength(true);

                if (!openNodes.Contains(otherNode)) // newly discovered node
                {
                    openNodes.Add(otherNode);
                }
                else if (tentativeGScore >= gScore[otherNode])
                {
                    continue; // This is not a better path
                }

                // Otherwise, this is the best path up til now
                cameFrom[otherNode] = currentNode;
                gScore[otherNode] = tentativeGScore;
                fScore[otherNode] = gScore[otherNode] + (endNode.GetXYZ() - otherNode.GetXYZ()).magnitude;
            }
        }

        // If endNode not reached, return empty list
        return shortestPath;
    }

    /// <summary>
    /// Determines whether there is a series of edges that continusously connect from the startNode to the endNode. 
    /// Optionally accepts a predicate function filtering out edges that do not meet the provided criteria.
    /// </summary>
    /// <param name="startNode"></param>
    /// <param name="endNode"></param>
    /// <param name="criteria"></param>
    /// <returns></returns>
    public bool DoesPathExist(Node startNode, Node endNode, NodeEdgePairs.EdgeFilterCriteria criteria = null)
    {
        HashSet<Node> closedSet = new HashSet<Node> { startNode };

        Queue<Node> queue = new Queue<Node>();
        queue.Enqueue(startNode);

        while (queue.Count > 0)
        {
            Node node = queue.Dequeue();
            foreach (Edge edge in edgesFromNode.GetEdges(node, criteria))
            {
                Node otherNode = edge.GetOtherNode(node);

                if (!closedSet.Contains(otherNode))
                {
                    if (criteria == null || criteria(edge))
                    {
                        if (otherNode == endNode)
                        {
                            return true;
                        }
                        queue.Enqueue(otherNode);
                        closedSet.Add(otherNode);
                    }
                }
            }
        }

        return false;
    }

    public struct RandomPath
    {
        public Node startNode;
        public Node finishNode;
        public List<Edge> pathEdges;
    }

    /// <summary>
    /// Generates a random path of the required length that does not intersect with itself. 
    /// Accepts an optional set of Nodes which should be ignored, or a dictionary of Edges to ignore 
    /// from a given node. Generates an error if required length cannot be met.
    /// </summary>
    /// <param name="requiredPathLength"></param>
    /// <param name="closedNodeSet"></param>
    /// <param name="closedEdgeSet"></param>
    /// <returns></returns>
    public RandomPath GetRandomPath(int requiredPathLength, HashSet<Node> closedNodeSet = null, Dictionary<Node, HashSet<Edge>> closedEdgeSet = null)
    {
        List<Edge> pathEdges = new List<Edge>();
        
        Node startNode = edgesFromNode.GetNodes().ToList()[Random.Range(0, nodes.Count)];

        closedNodeSet = closedNodeSet == null ? new HashSet<Node> { startNode } : closedNodeSet;
        closedEdgeSet = closedEdgeSet == null ? new Dictionary<Node, HashSet<Edge>>() : closedEdgeSet;

        Node currentNode = startNode;
        while (pathEdges.Count < requiredPathLength)
        {
            if (!closedEdgeSet.Keys.Contains(currentNode))
            {
                closedEdgeSet.Add(currentNode, new HashSet<Edge>());
            }
            HashSet<Edge> availableEdges = new HashSet<Edge>(
                            edgesFromNode.GetEdges(currentNode).Where(x => !closedEdgeSet[currentNode].Contains(x)
                                                                        && !closedNodeSet.Contains(x.GetOtherNode(currentNode))
                                                                        && !pathEdges.Contains(x)).ToList());
            // Back-track if we've reached a dead end
            if (availableEdges.Count == 0)
            {
                if (currentNode == startNode)
                {
                    Debug.LogError("Random Path could not be found for start node: " + startNode + " with target path length: " + requiredPathLength);
                }
                Edge previousPath = pathEdges[pathEdges.Count - 1];
                Node previousNode = previousPath.GetOtherNode(currentNode);
                closedEdgeSet[previousNode].Add(previousPath);
                pathEdges.Remove(previousPath);
                currentNode = previousNode;
                continue;
            }
            Edge nextEdge = availableEdges.ToList<Edge>()[Random.Range(0, availableEdges.Count)];
            pathEdges.Add(nextEdge);
            closedNodeSet.Add(currentNode);
            currentNode = nextEdge.GetOtherNode(currentNode);
        }

        RandomPath randomPath = new RandomPath();
        randomPath.startNode = startNode;
        randomPath.finishNode = currentNode;
        randomPath.pathEdges = pathEdges;

        return randomPath;
    }
}

public class NodeEdgePairs
{
    Dictionary<Node, List<Edge>> nodeEdgePairs;

    public NodeEdgePairs()
    {
        this.nodeEdgePairs = new Dictionary<Node, List<Edge>>();
    }

    public void Add(Node node, Edge edge)
    {
        if (!nodeEdgePairs.Keys.Contains(node))
        {
            nodeEdgePairs.Add(node, new List<Edge> { edge });
        }
        else if (!nodeEdgePairs[node].Contains(edge))
        {
            nodeEdgePairs[node].Add(edge);
        }
    }

    public HashSet<Node> GetNodes()
    {
        return new HashSet<Node>(nodeEdgePairs.Keys);
    }

    public delegate bool EdgeFilterCriteria(Edge edge);

    /// <summary>
    /// Returns a list of the edges emanating from the given node. Optionally accepts  
    /// a predicate function filtering out edges that do not meet the provided criteria.
    /// </summary>
    /// <param name="node"></param>
    /// <param name="criteria"></param>
    /// <returns></returns>
    public List<Edge> GetEdges(Node node, EdgeFilterCriteria criteria = null)
    {
        List<Edge> edges = nodeEdgePairs[node];
        if (criteria == null)
        {
            return edges;
        }
        else
        {
            return edges.Where(x => criteria(x)).ToList();
        }
    }

    public delegate int EdgeComparison(Edge edge1, Edge edge2);

    /// <summary>
    /// Sorts the list of edges for the provided node given a comparison function.
    /// </summary>
    /// <param name="node"></param>
    /// <param name="compare"></param>
    public void SortEdges(Node node, EdgeComparison compare)
    {
        nodeEdgePairs[node].Sort((e1, e2) => compare(e1, e2));
    }

    /// <summary>
    /// Sorts all edges maintained in the object based on the given comparison function.
    /// </summary>
    /// <param name="compare"></param>
    public void SortAllEdges(EdgeComparison compare)
    {
        foreach (KeyValuePair<Node, List<Edge>> pair in nodeEdgePairs)
        {
            pair.Value.Sort((e1, e2) => compare(e1, e2));
        }
    }
}

public class Edge
{
    Node node1;
    Node node2;
    List<Node> subNodes;
    float weight;

    public Edge(Node node1, Node node2, List<Node> subNodes = null, float weight = 1f)
    {
        this.node1 = node1;
        this.node2 = node2;
        this.weight = weight;

        // Ensure that no duplicates are added
        // Note: documentation states that original order is not guaranteed,
        //      however, implementation and testing show that it is currently preserved.
        //      (I.e. preservation of order may be subject to change with future versions of .NET)
        this.subNodes = subNodes.Distinct().ToList();
    }

    public Node GetFirstNode()
    {
        return node1;
    }

    public Node GetLastNode()
    {
        return node2;
    }

    /// <summary>
    /// Returns a list of the primary start and end nodes in order. (Does not include sub-nodes)
    /// </summary>
    /// <returns></returns>
    public List<Node> GetNodes()
    {
        return new List<Node> { node1, node2 };
    }

    /// <summary>
    /// Returns the other node provided one of the nodes. (Ignores sub-nodes)
    /// </summary>
    /// <param name="node"></param>
    /// <returns></returns>
    public Node GetOtherNode(Node node)
    {
        return node == node1 ? node2 : node1;
    }

    /// <summary>
    /// Returns the other node provided one of the nodes. (Excludes primary start and end nodes)
    /// </summary>
    /// <param name="node"></param>
    /// <returns></returns>
    public List<Node> GetSubNodes()
    {
        return subNodes;
    }

    /// <summary>
    /// Sets the weight value for this Edge. For best results, use values >= 1.
    /// </summary>
    /// <param name="weight"></param>
    public void SetWeight(float weight)
    {
        this.weight = weight;
    }

    /// <summary>
    /// Returns the total distance between each Node in the Edge (including sub-nodes).
    /// If optional weighted parameter is set to true, will multiply distance by weight.
    /// </summary>
    /// <param name="weighted"></param>
    /// <returns></returns>
    public float GetPathLength(bool weighted = false)
    {
        float length = 0;
        List<Vector3> path = this.GetPath();
        for (int i = 0; i < path.Count - 1; i++)
        {
            length += Vector3.Distance(path[i], path[i + 1]);
        }
        return weighted ? weight * length : length;
    }

    /// <summary>
    /// Returns the total distance between each position in the provided path.
    /// If optional weighted parameter is set to true, will multiply distance by weight.
    /// </summary>
    /// <param name="path"></param>
    /// <param name="weighted"></param>
    /// <returns></returns>
    public float GetPathLength(List<Vector3> path, bool weighted = false)
    {
        float length = 0;
        for (int i = 0; i < path.Count - 1; i++)
        {
            length += Vector3.Distance(path[i], path[i + 1]);
        }
        return weighted ? weight * length : length;
    }

    /// <summary>
    /// Returns a list of the Vector3 coordinates for each node in the edge, 
    /// including any sub-nodes, 
    /// starting at the provided startNode.
    /// If startNode is not provided, will start from first node of Edge
    /// </summary>
    /// <param name="startNode"></param>
    /// <returns></returns>
    public List<Vector3> GetPath(Node startNode)
    {
        if (startNode != node1 && startNode != node2)
        {
            Debug.LogError("The provided startNode (" + startNode + ") is not one of the primary Nodes of this Edge.");
        }
        List<Vector3> path = new List<Vector3>();
        path.Add(startNode == node1 ? node1.GetXYZ() : node2.GetXYZ());
        if (subNodes != null)
        {
            foreach (Node subNode in startNode == node1 ? subNodes : subNodes.AsEnumerable().Reverse())
            {
                path.Add(subNode.GetXYZ());
            }
        }
        path.Add(startNode == node1 ? node2.GetXYZ() : node1.GetXYZ());
        return path;
    }

    /// <summary>
    /// Returns a list of the Vector3 coordinates for each node in the edge, including any sub-nodes.
    /// </summary>
    /// <returns></returns>
    public List<Vector3> GetPath()
    {
        List<Vector3> path = new List<Vector3>();
        path.Add(node1.GetXYZ());
        if (subNodes != null)
        {
            foreach (Node subNode in subNodes)
            {
                path.Add(subNode.GetXYZ());
            }
        }
        path.Add(node2.GetXYZ());
        return path;
    }

    /// <summary>
    /// Provides an interpolated point from the startNode to the other Node 
    /// given t, where the startNode is at t = 0.0 and the other node is at t = 1.0. 
    /// If no startNode provided, defaults to the first node of this Edge.
    /// </summary>
    /// <param name="startNode"></param>
    /// <param name="t"></param>
    /// <returns></returns>
    public Vector3 GetPointAlongEdge(float t, Node startNode = null)
    {
        startNode = startNode == null ? node1 : startNode;
        if (startNode != node1 && startNode != node2)
        {
            Debug.LogError("The provided startNode is not one of the primary Nodes of this Edge.");
        }

        List<Vector3> path = GetPath(startNode);
        float totalLength = GetPathLength(path);

        float accumulatedT = 0f;
        float accumulatedDistance = 0f;
        for (int i = 0; i < path.Count - 1; i++)
        {
            Vector3 currentVector = path[i + 1] - path[i];
            float currentVectorLength = currentVector.magnitude;
            if (accumulatedDistance + currentVectorLength > totalLength)
            {
                return path[i] + (t - accumulatedT) * currentVector;
            }
            else
            {
                accumulatedDistance += currentVectorLength;
                accumulatedT += currentVectorLength / t;
            }
        }
        // If a point has not yet been found, return the last point on the path
        return path[path.Count - 1];
    }

    /// <summary>
    /// Returns a new Edge with all nodes (including any sub-nodes) in the reverse order.
    /// </summary>
    /// <returns></returns>
    public Edge GetReversed()
    {
        List<Node> reversedSubNodes = null;
        if (subNodes != null)
        {
            reversedSubNodes = subNodes.AsEnumerable().Reverse().ToList();
        }
        return new Edge(node2, node1, reversedSubNodes);
    }
}

public class Node
{
    Vector3 coordinates;

    public Node(Vector3 coordinates)
    {
        this.coordinates = coordinates;
    }

    public Node(Vector2 coordinates)
    {
        this.coordinates = new Vector3(coordinates.x, coordinates.y, 0f);
    }

    public Vector3 SetCoordinates(Vector3 coords)
    {
        coordinates = coords;
        return coordinates;
    }

    public Vector3 GetXYZ()
    {
        return coordinates;
    }

    public Vector2 GetXY()
    {
        return coordinates;
    }

    public Vector2 GetXZ()
    {
        return new Vector2(coordinates.x, coordinates.z);
    }
}
