using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class Graph
{
    HashSet<Edge> edges = new HashSet<Edge>();
    NodeToEdges nodeToEdges = new NodeToEdges();

    public void AddUndirectedEdge(Edge newEdge)
    {
        AddDirectedEdge(newEdge);
        AddDirectedEdge(newEdge.GetReversed());
    }

    public void AddDirectedEdge(Edge newEdge)
    {
        edges.Add(newEdge);

        foreach (Node node in newEdge.GetNodes())
        {
            nodeToEdges.Add(node, newEdge);
        }
    }

    public delegate bool EdgeFilterCriteria(Edge edge);

    public HashSet<Edge> GetEdges(Node node, EdgeFilterCriteria criteria = null)
    {
        HashSet<Edge> edges = nodeToEdges.GetEdges(node);
        if (criteria == null)
        {
            return edges;
        }
        else
        {
            return new HashSet<Edge>(edges.ToList().Where(x => criteria(x)));
        }
    }
}

public class NodeToEdges
{
    Dictionary<Node, HashSet<Edge>> nodeToEdges;

    public NodeToEdges()
    {
        this.nodeToEdges = new Dictionary<Node, HashSet<Edge>>();
    }

    public void Add(Node node, Edge edge)
    {
        if (!nodeToEdges.Keys.Contains(node))
        {
            nodeToEdges.Add(node, new HashSet<Edge> { edge });
        }
        else
        {
            nodeToEdges[node].Add(edge);
        }
    }

    public HashSet<Edge> GetEdges(Node node)
    {
        return nodeToEdges[node];
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
    /// Sets the weight value for this Edge.
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
            Debug.LogError("The provided startNode is not one of the primary Nodes of this Edge.");
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
