import string, sys, random, uuid

class Graph:
  # Use directed graph by default
  def __init__(self,directionless=False):
    self.__edges = {}
    self.__nodes = {}
    self.__directionless=directionless
    self.__parent_to_child = {} #keyed by node
    self.__child_to_parent = {}

  def __str__(self):
    return self.get_report()

  def get_report(self):
    ostr = ''
    ostr += "Nodes: "+str(len(self.__nodes.keys()))+"\n"
    ostr += "Edges: "+str(len(self.__edges.keys()))+"\n"
    return ostr
    
  def get_edges(self):
    return self.__edges.values()

  # get a edges given a node
  def get_node_edges(self,node,type="both"):
    if type == "both":
      return [self.__edges[x] for x in self.__edges if node.get_id() in self.__edges[x].get_node_ids()]
    elif type == "outgoing":
      return [self.__edges[x] for x in self.__edges if node.get_id() == self.__edges[x].get_node_ids()[0]]
    elif type == "incoming":
      return [self.__edges[x] for x in self.__edges if node.get_id() == self.__edges[x].get_node_ids()[1]]
    sys.stderr.write("ERROR: type is not a type "+type+"\n")
    sys.exit()

  def get_nodes(self):
    return self.__nodes.values()
    #return self.__nodes

  def add_node(self,node):
    #if node.has_edges():
    #  sys.stderr.write("ERROR: nodes can only be added to a graph before edges are set\n")
    #  sys.exit()
    self.__nodes[node.get_id()] = node
    return

  def get_children(self,node):
    if self.find_cycle() or self.__directionless:
      sys.stderr.write("ERROR: do cannot find a branch when there are cycles in the graph\n")
      sys.exit()
    v = self.__get_children(node.get_id())
    return [self.__nodes[i] for i in v]

  #assumes you have no cycles and its a directional graph
  def __get_children(self,nodeid):
    if nodeid not in self.__parent_to_child: #this is a leaf so traverse no farther
      return []
    kids = []
    for j in self.__parent_to_child[nodeid]:
      v = self.__get_children(j)
      for k in v:  
         if k not in kids: kids.append(k)
      kids.insert(0,j)
    return kids

  def get_roots(self):
    if self.__directionless:
      sys.stderr.write("ERROR: can't get roots of an undirected graph\n")
      sys.exit()
    outputids = self.__nodes.keys()
    #print outputids
    rootset  = set(outputids) -  set(self.__child_to_parent.keys())
    return [self.__nodes[x] for x in rootset]

  def add_edge(self,edge,verbose=True):
    #make sure nodes are in the nodes
    if edge.get_node1().get_id() not in self.__nodes:
      sys.stderr.write("ERROR: node should be in graph\n")
      sys.exit()
      #self.__nodes[edge.get_node1().get_id()] = edge.get_node1()
    if edge.get_node2().get_id() not in self.__nodes:
      sys.stderr.write("ERROR: node should be in graph\n")
      sys.exit()
      #self.__nodes[edge.get_node2().get_id()] = edge.get_node2()
    # now add edge
    id = edge.get_id()
    #sys.stderr.write(id+"\n")
    if id in self.__edges:
      sys.stderr.write("WARNING edge is already there. not adding again\n")
      return
    self.__edges[id] = edge
    ids = edge.get_node_ids()    
    if ids[0] not in self.__parent_to_child:
      self.__parent_to_child[ids[0]] = {}
    if ids[1] in self.__parent_to_child[ids[0]] and verbose==True:
      if verbose: sys.stderr.write("Warning repeat edge.\n")
      return
    self.__parent_to_child[ids[0]][ids[1]] = edge.get_id()
    #if edge.is_directionless():
    if ids[1] not in self.__child_to_parent:
      self.__child_to_parent[ids[1]] = {}
    if ids[0] in self.__child_to_parent[ids[1]] and verbose == True:
      sys.stderr.write("WARNING overwriting repeat edge.\n")
    self.__child_to_parent[ids[1]][ids[0]] = edge.get_id()
    return

  def get_status_string(self):
    ostr = ''
    ostr += "----------------\n"
    ostr += "Node count: "+str(len(self.__nodes.keys()))+"\n"
    ostr += "Edge count: "+str(len(self.__edges.keys()))+"\n"
    if not self.__directionless:
      ostr += "Root count: "+str(len(self.get_roots()))+"\n"
    return ostr

  def remove_node(self,node):
    nid = node.get_id()
    #remove edges associated with this node
    edges = self.get_node_edges(node,type="both")
    for e in edges:
      self.remove_edge(e)
    del self.__nodes[nid]

  # remove edge
  def remove_edge(self,edge):
    #sys.stderr.write(str(len(self.__edges.keys()))+" edges\n")
    #sys.stderr.write(str(self.__edges.keys())+" edges\n")
    #sys.stderr.write(str(self.get_report())+" remove edge\n")
    #sys.stderr.write(str([x.get_node_ids() for x in self.__edges.values()])+" remove edge\n")
    if edge.get_id() not in self.__edges:
      sys.stderr.write("WARNING: edge already removed\n")
      return
    nodeids = edge.get_node_ids()
    node1 = self.__nodes[nodeids[0]]
    node2 = self.__nodes[nodeids[1]]
    edges_to_remove = set()
    if node1.get_id() in self.__parent_to_child:
      if node2.get_id() in self.__parent_to_child[node1.get_id()]:
        edges_to_remove.add(self.__parent_to_child[node1.get_id()][node2.get_id()])
        del self.__parent_to_child[node1.get_id()][node2.get_id()]
      if len(self.__parent_to_child[node1.get_id()]) == 0:
        del self.__parent_to_child[node1.get_id()]
    if node2.get_id() in self.__child_to_parent: 
        if node1.get_id() in self.__child_to_parent[node2.get_id()]:
          edges_to_remove.add(self.__child_to_parent[node2.get_id()][node1.get_id()])
          del self.__child_to_parent[node2.get_id()][node1.get_id()]
        if len(self.__child_to_parent[node2.get_id()]) == 0:
          del self.__child_to_parent[node2.get_id()]
    if len(edges_to_remove) == 0:
      sys.stderr.write("WARNING no edges removed\n")
    #sys.stderr.write(str(len(self.__edges.keys()))+" edges\n")
    #sys.stderr.write(str(self.__edges.keys())+" edges\n")
    #sys.stderr.write(str(edges_to_remove)+" edges\n")
    for eid in edges_to_remove:
      del self.__edges[eid]

  #remove cycles by mergine cyclic nodes into single nodes
  #their payloads are added to a list
  def merge_cycles(self):
    #delete any self cycles first
    for i in self.__parent_to_child:
      for j in self.__parent_to_child[i]:
        if self.__parent_to_child[i]==j: 
          self.remove_edge(self.__edges[self.__parent_to_child[i][j]])
    while True:
      res = self.find_cycle()
      if not res:  return # we have finished.. there are no cycles
      # If we are here we need to merge
      if len(res) == 1: 
        sys.stderr.write("ERROR: Unexpected Self-cycle.\n")
        sys.exit()
      resids = [x.get_id() for x in res]
      # merge edges unless that edge is to one of the nodes we are removing
      for i in range(1,len(res)):
        for v in res[i].get_payload(): res[0].get_payload().append(v)
        if res[i].get_id() in self.__parent_to_child:
          for e2id in self.__parent_to_child[res[i].get_id()]:
            if e2id not in resids:
              nedge = Edge(res[0],res[i])
              #sys.stderr.write("adding :"+nedge.get_id()+"\n")
              if res[0].get_id() not in self.__parent_to_child:
                self.add_edge(nedge,verbose=False)
              elif res[1].get_id() not in self.__parent_to_child[res[0].get_id()]:
                self.add_edge(nedge,verbose=False)
        #      #if self.__directionless:
        #      #  sys.stderr.write('adding_edge2'+"\n")
        #      #  self.add_edge(Edge(res[i],res[0]),verbose=False)
      # remove any nodes and edges connected to nodes we are removing
      for r in res[1:]:
        self.remove_node(r)

  #return a single cycle, greedy first one found
  #in terms of nodes return as an array of nodes or None
  def find_cycle(self):
    #depth first search through nodes
    for nid in self.__nodes.keys():
      res = self.__find_cycle_node([],nid)
      if res: return [self.__nodes[x] for x in res]
    return None

  # From some node
  def get_directed_paths_from_node(self,node,prev=[]):
    if self.__directionless:
      sys.stderr.write("ERROR: Can't find paths from directionless graph\n")
      sys.exit()
    if self.find_cycle(): 
      sys.stderr.write("ERROR: Can't find paths when a cycle is present.\n")
      sys.exit()
    id = node.get_id()
    nprev = prev[:]
    nprev.append(node)
    if id in self.__parent_to_child:
      output = []
      for nextnode_id in self.__parent_to_child[id]:
        vs = self.get_directed_paths_from_node(self.__nodes[nextnode_id],nprev)
        if vs: 
          for v in vs: output.append(v)
      return output
    else:
      # we are at a leaf
      return nprev

  #Internal function
  # Return the first cycle found form a starting node in terms of an
  # array of node ids          
  def __find_cycle_node(self,starting_ids,current):
    if current in starting_ids: return starting_ids
    if current not in self.__parent_to_child: return None
    for i in self.__parent_to_child[current]:
      newstarts = starting_ids[:]
      newstarts.append(current)
      res = self.__find_cycle_node(newstarts,i)
      if res: return res
    return None

  def partition_graph(self,verbose=False):
    visited_nodes = set()
    ns = self.get_nodes()
    #exclude_ids = exclude_ids | set([x.get_id() for x in g.get_nodes()])
    node_sets = []
    z = 0
    for n in ns:
      z += 1
      if verbose: sys.stderr.write("partitioning: "+str(z)+'/'+str(len(ns))+"       \r")
      nids = set([y.get_id() for y in self.connected_nodes(n,exclude_ids=visited_nodes)])
      node_sets.append(nids)
      visted_nodes = visited_nodes | nids
    if verbose: sys.stderr.write("\n")
    #node_sets = [set([y.get_id() for y in self.connected_nodes(x)]) for x in ns]
    results = []
    tot = len(node_sets)
    while len(node_sets) > 0:
      if verbose: sys.stderr.write('reducing node sets '+str(tot-len(node_sets)+1)+'/'+str(tot)+"       \r")
      v = node_sets.pop()
      #print 'v'
      #print v
      added = False
      for i in range(0,len(results)):
        if v & results[i]:
          #print results[i]
          results[i] = results[i] | v
          added = True
          break
      if not added: results.append(v)
    if verbose: sys.stderr.write("\n")
    g_results = []
    z = 0
    for r in results:
      z += 1
      if verbose: sys.stderr.write("making graph: "+str(z)+"/"+str(len(results))+"       \r")
      g = Graph(directionless=self.__directionless)
      for nid in r:
        g.add_node(self.__nodes[nid])
        for e in self.get_node_edges(self.__nodes[nid],type="outgoing"):
          g.add_node(e.get_node2())
          g.add_edge(Edge(e.get_node1(),e.get_node2()),verbose=False)
      #for v in [x.get_payload() for x in g.get_nodes()]:
      #  for txs in v:
      #es    print txs.get_gene_name()
      g_results.append(g)
    return g_results

  def connected_nodes(self,node,exclude_ids=None):
    r =  _depth_traverse_node(self,node,visited=exclude_ids)
    if not r: return []
    #sys.stderr.write(str(node.get_id())+"\n")
    #sys.stderr.write(str(self.__nodes.keys())+"\n")
    #sys.stderr.write("\n"+self.get_report()+"\n")
    #sys.stderr.write(str(self.__child_to_parent.keys())+"\n")
    #sys.stderr.write(str([x for x in r])+"\n")
    return [self.__nodes[x] for x in r]

def _depth_traverse_node(g,node,visited=None):
  #print visited
  if not visited: visited = set()
  if node.get_id() in visited:
    return visited
  visited.add(node.get_id())
  #ns = [x for x in g.get_nodes() if x.get_id() not in visited]
  #if len(ns) == 0: return
  es = g.get_node_edges(node)
  if not es: return visited
  tot = set()
  tot = tot | visited
  for e in es:
    n2 = e.get_node2()
    if n2.get_id() in tot: continue
    v = _depth_traverse_node(g,n2,tot)
    tot = tot | v
  return tot

# directed graph by default
class Edge:
  def __init__(self,node1,node2,directionless=False,weight=None):
    self.__node1 = node1
    self.__node2 = node2
    self.__directionless = directionless
    self.__weight = weight
    self.__id = str(uuid.uuid4())
  def set_weight(self,weight):  self.__weight = weight
  def get_weight(self): return self.__weight
  def get_node_ids(self): return [self.__node1.get_id(),self.__node2.get_id()]
  def is_directionless(self): return self.__directionless
  def get_node1(self): return self.__node1
  def get_node2(self): return self.__node2
  def get_id(self): return self.__id

#payload is a list. When nodes get merged lists are concatonated.
class Node:
  def __init__(self,payload=None):
    self.__id = str(uuid.uuid4())
    self.__payload = []
    #self.__incoming_edges = {}
    #self.__outgoing_edges = {}
    if payload != None:
      self.__payload = payload
  #def has_edges(self):
  #  if len(self.__incoming_edges.keys()) > 0: return True
  #  if len(self.__outgoing_edges.keys()) > 0: return True
  #  return False
  def get_payload(self):
    return self.__payload
  def set_payload(self,payload):
    self.__payload = payload
  def get_id(self):
    return self.__id
  #def set_child_node(self,node2):
  #  e = Edge(self,node2)
  #  self.__outgoing_edges[e.get_id()] = e
  #def set_parent_node(self,node2):
  #  e = Edge(node2,self)
  #  self.__incoming_edges[e.get_id()] = e
