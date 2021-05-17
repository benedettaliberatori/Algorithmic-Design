from numbers import Number
from collections import defaultdict
import numpy as np
from binheap_dijkstra import binheap


inf = float('inf')

def pair_order(a,b) -> bool:
    '''
    Compares two pairs w.r.t.
    the second element of both
    and the '<=' operation. 

    '''
    a1, a2 = a
    b2, b2 = b
    return a2 <= b2

'''-----------------------------------------------------------------------------------------------------------------'''
class graph:  
    '''A graph class

    Members
    -------
    adj: 
        A dictionary of dictionaries. 
        -The keys of the outer dictionary are the nodes, 
        -the values are dictionaries with: 
                -keys that are the neighbours nodes 
                -values that are the weights associated.
    
    
    ----------
    undirected: 
        A boolean stating whether 
        the graph is undirected or not.
        Default is False. 

    
    '''

    def __init__(self, undirected: bool=False):
        self.adj = defaultdict(dict)
        self.undirected = undirected
    
    @property 
    def nodes(self):
        '''
        Returns a set containing
        the nodes in the graph.
        '''
        return {k for k in self.adj}
    

    def insert_edge(self: int, src: int, dst: int , weight: int = None) -> None:

        ''' 
        Inserts an edge into the graph.
        Inserts an additional edge from 
        dst to src with the same weight 
        if the graph is undirected. 
            
        Parameters
        ----------
        src: int
            The source node id.
        
        dst: int
            The destination node id. 
        
        weight: int 
            The weight of the edge. 
            Default is None and an edge 
            with weight 1 is inserted. 
       
        '''
        self.adj[src][dst] = weight if weight else 1
        if self.undirected:
            self.adj[dst][src] = weight if weight else 1 
        
        if dst not in self.adj:
            self.adj[dst]={}
        
    
    def print_graph(self):
        '''
        Example : 
        1: {6: 1, 5: 1}
        6: {8: 1}
        5: {1: 3, 6: 1}
        8: {1: 1, 7: 1}
        7: {8: 1}
        4: {8: 3, 7: 1}
        3: {4: 3, 2: 1}
        2: {3: 2}
        '''
        for item in self.adj:
            print("{}: {}".format(item,self.adj[item]))


'''--------------------------------------------------------------------------------------'''

def relax(H: binheap, u: int, v: int, w: int, dist: list, pred: list) -> None:
    ''' 
        Updates the distance of the node v 
        if needed.
                            
        Parameters
        ----------
        H: binheap
            Binary Heap implementing
            the Priority Queue. 
        
        u: int
           node id
            
        v: int
           node id

        w: int 
           weight of the edge
           from u to v
           
        dist: list
            Stores in position i-1 
            the distance of the node
            i from src. 
            
    
        pred: list
            Stores in position i-1 
            the predecessor of the node
            i in the shortest path

                
        '''

    if dist[u-1] + w < dist[v-1]:
        
        dist[v-1]= dist[u-1] + w
        H.decrease_key(v, (v,dist[v-1])) 
        pred[v-1]=u
    


def Dijkstra(G: graph, src: int):
    ''' 
        Implements the Dijkstra's algorithm for single
        source shortest paths, using binary heaps.

                   
        Parameters
        ----------
        G: graph
        
        src: int
            Id of the node selected
            as the source of the 
            Dijkstra algorithm. 
        
        
        
        Returns 
        ----------
        
        dist:
        a list where the (i-1)-th
        element is the shortest distance 
        of the i-th node from src

        pred: 
        a list where the (i-1)-th
        element is the predecessor of 
        the i-th node in the shortest 
        path from src.
        
        
        
    '''

    assert src in G.nodes , 'Source node not found'

    # Initialization
    dist=[inf for i in G.nodes]
    pred=[None for i in G.nodes]
    dist[src-1] = 0

    # Filling the Priority Queue with couples
    # of node and distance from the src (up to now).
    # using the pair_order operation for 
    # comparison.  
    H = binheap([(v, dist[v-1]) for v in G.nodes], pair_order)

    while not H.is_empty():
        u , _ = H.remove_minimum()
        
        
        
        for v in G.adj[u]:
               
                relax(H, u, v, G.adj[u][v], dist, pred)
                
    
    return dist, pred 


def witness_search(G: graph, src: int, dst: int, maxdist: int ,x=None) -> int:
    
    assert src in G.nodes , 'Source node not found'

    

    # Initialization
    dist =[inf for n in G.nodes if n!=x]
    pred = [None for n in G.nodes if n!=x]
    dist[src-1] = 0

    H = binheap([(v, dist[v-1]) for v in G.nodes], pair_order)
    

    while not H.is_empty():
        u , u_dist = H.remove_minimum()
       
        if u == dst: return u_dist
         
        for v in G.adj[u]:
            if v not in x:
                if G.adj[u][v] > maxdist: break
                                        
                else:
                    relax(H, u, v, G.adj[u][v], dist, pred)
            


     

def shortcut(G: graph, v: int, removed: set)-> None:
    ''' 
        Adds to the graph G the shortcuts 
        needed for the contraction
        of node v.

                   
        Parameters
        ----------
        G: graph

        v: int
           node to be contracted.

        removed: set 
           contains the nodes
           already contracted.
        
    '''

    assert v in G.nodes , 'Node not found'

    
    # Retrieve the predecessors of v
    # among those not already contracted 
    pred = set()
    for p in G.nodes:
        for k in G.adj[p]:
            if k ==  v and p not in removed:
                pred.add(p)
    
    # For each couple of predecessor-successor of v
    for u in pred:
        for w in G.adj[v]:
            if w not in removed:
                shortcut_cost = G.adj[u][v] + G.adj[v][w]
                
                # compute the shortest distance between u and w 
                witness_cost = witness_search(G, u, w, shortcut_cost,removed)
            
            
            # add shortcut if needed
                if shortcut_cost < witness_cost:
                    G.insert_edge(u,w,shortcut_cost)

def preprocessing(G: graph)->None:
    '''
    Decorates the graph G 
    with all the shortcuts.

    '''
    removed=set()
    for x in G.nodes:
        removed.add(x)
        shortcut(G,x, removed)
    

def Bid_Dijkstra(G: graph, src: int, dst: int)-> int:
    '''
    Performs a bidirectional version of Dijkstra
    on an augmented graph to compute the shortest
    path distance between src and dst.

    Parameters
        ----------
        G: graph

        src: int
            Source node
        
        dst: int
            Destination node
    
    Returns 
        ----------

        Computed minimum distance
    

    '''
    
    #Upward graph
    U = graph()
    
    #Downward graph 
    D = graph()

    # Filling both with all the nodes in G
    for n in G.nodes: 
        U.adj[n]={}
        D.adj[n]={}
      
    for i in G.adj:
        for j in G.adj[i]: 
            if j > i: 
                # insert edges from less to more important
                U.insert_edge(i,j,G.adj[i][j])
            else:
                # insert edges from more to less important,
                # inverted ( to perform backward Dijkstra)
                D.insert_edge(j,i,G.adj[i][j]) 
    
        
    
    # Initialization
    df = [inf for n in U.nodes]
    db = [inf for n in D.nodes]
    df[src-1] = 0
    db[dst-1] = 0
    pf = [None for n in U.nodes]
    pb = [None for n in D.nodes]

    
    # Filling two distinct Priority Queues
    Qf = binheap([(v, df[v-1]) for v in U.nodes], pair_order)
    Qb = binheap([(v, db[v-1]) for v in D.nodes], pair_order)
    Dist = inf


    while not Qf.is_empty() and not Qb.is_empty():
        u , _ = Qf.remove_minimum()
        v , _ = Qb.remove_minimum()
        
                

        for i in U.adj[u]:
            relax(Qf, u, i, U.adj[u][i], df, pf)
            
            if db[i-1] is not inf:
                Dist = min(Dist, df[u-1] + U.adj[u][i] + db[i-1] )
                
        
        for j in D.adj[v]:
            relax(Qb, v, j, D.adj[v][j], db, pb)
            
            
            if df[i-1] is not inf:
                Dist = min(Dist, db[v-1] + D.adj[v][j] + df[j-1])
                
        
    return Dist
    
    
    
if __name__ == '__main__':
    g = graph()
    
    g.insert_edge(1, 6, 1)
    g.insert_edge(5, 1, 3)
    g.insert_edge(1, 5, 1)
    g.insert_edge(5, 6, 1)
    g.insert_edge(6, 8, 1)
    g.insert_edge(8, 1, 1)
    g.insert_edge(8,7,1)
    g.insert_edge(7,8,1)
    g.insert_edge(4,8,3)
    g.insert_edge(4,7,1)
    g.insert_edge(3,4,3)
    g.insert_edge(3,2,1)
    g.insert_edge(2,3,2)
    
    print("Input graph:\n")
    g.print_graph()
    print("\n")
    
    src = 2
    print(f"Result of Dijkstra(g, {src}): { Dijkstra(g, src)}\n")

    preprocessing(g)
    print("Augmented graph:\n")
    g.print_graph()
    print("\n")
    dst = 6
    print(f"Result of Bid_Dijkstra(g, {src}, {dst}): {Bid_Dijkstra(g, src, dst)}\n")
   

    
    