#include <fstream>
#include <iostream>
#include <vector>
#include <bits/stdc++.h>
using namespace std;
#define MAXV 10000 
#define oo 0x3f3f3f3f 
ofstream outputStream;
int previo[ MAXV ];
void reverse(vector<int> &path, int start, int end, int n)
{
	while (end - start > 0)
	{
		
		int temp = path[start % n];
		path[start % n] = path[end % n];
		path[end % n] = temp;
		start++;
		end--;
	}
};

int is_path_shorter(long int **graph, int v1, int v2, int v3, int v4, int &total_dist)
{
	if ((graph[v1][v3] + graph[v2][v4]) < (graph[v1][v2] + graph[v3][v4]))
	{
		total_dist -= (graph[v1][v2] + graph[v3][v4] - graph[v1][v3]
				- graph[v2][v4]);
		return 1;
	}
	return 0;
};

int twoOpt(long int **graph, vector<int> &path, int &pathLength, int n)
{
	int counter = 0;
	int term_cond = 5;
	int old_distance = pathLength;
	
	int v1, v2, u1, u2;

	
	for (int i = 0; i < n; i++)
	{
		
		u1 = i;
		v1 = (i+1) % n; 

		
		for (int j = i + 2; (j + 1) % n != i; j++)
		{
			
			u2 = j % n;
			v2 = (j+1) % n;

		
			if (is_path_shorter(graph, path[u1], path[v1], path[u2],
					path[v2], pathLength))
			{

				
				reverse(path, i + 1, j, n); // v1, u2

				if (pathLength - old_distance < 5 && counter == term_cond)
					break;

				
				if (pathLength - old_distance > term_cond )
					i = 0;

				old_distance = pathLength;
				counter++;
			}
		}
	}
	return pathLength;
};

int get_path_length(int **graph, vector<int> &path, int size){
	int length = 0;
	for (int i = 0; i < size-1; i++)
	{
		length += graph[path[i]][path[i+1]];
	}
	length += graph[path[size-1]] [path[0]]; 
	return length;
};
struct Arista{
	int nodo; 
	int costo; 
	Arista(int _nodo, int _costo) : nodo(_nodo), costo(_costo) {} 
	Arista() : nodo(-1), costo(-1) {} 
};

struct Grafo{
	vector<Arista> G[MAXV];
	int V; 
	int A;
};

struct Nodo{
	int nodo; 
	int costo; 
	Nodo(int _nodo, int _costo) : nodo(_nodo), costo(_costo) {} 
	bool operator <(const Nodo &b) const {
		return costo > b.costo;
	}
};
int a;
void camino( int destino , vector<int> previo ){ 
	if(previo[destino]!= -1)
		camino(previo[destino],previo);
	if(a!=destino)
		outputStream<<destino<<endl;	
}

int dijkstraCOSTO(const int inicial, const int destino, const Grafo grafo){
	priority_queue<Nodo> Q;
	vector<int> distancia(grafo.V, oo);
	vector<bool> visitado(grafo.V, false);
	distancia[inicial] = 0;
	Q.push(Nodo(inicial, 0));
	while(!Q.empty()) {
		Nodo actual = Q.top(); 
        Q.pop();
		visitado[actual.nodo] = true;
		if (actual.nodo == destino)
			return actual.costo;
		int T = (int)grafo.G[actual.nodo].size();
		for(int i = 0; i < T; ++i){
			if (!visitado[grafo.G[actual.nodo][i].nodo] && ((distancia[actual.nodo] + grafo.G[actual.nodo][i].costo) < distancia[grafo.G[actual.nodo][i].nodo])){
				distancia[grafo.G[actual.nodo][i].nodo] = actual.costo + grafo.G[actual.nodo][i].costo;
				Q.push(Nodo(grafo.G[actual.nodo][i].nodo, actual.costo + grafo.G[actual.nodo][i].costo));
			}
		}
	}
	return -1;
}

int dijkstraCAMINO(const int inicial, const int destino, const Grafo grafo){
	priority_queue<Nodo> Q;
	vector<int> distancia(grafo.V, oo);
	vector<bool> visitado(grafo.V, false); 
    vector<int> previo(grafo.V, -1);
	distancia[inicial] = 0;
	Q.push(Nodo(inicial, 0));
	while(!Q.empty()) {
		Nodo actual = Q.top(); 
        Q.pop();
		visitado[actual.nodo] = true;
		int T = (int)grafo.G[actual.nodo].size();
		for(int i = 0; i < T; ++i){
			if (!visitado[grafo.G[actual.nodo][i].nodo] && ((distancia[actual.nodo] + grafo.G[actual.nodo][i].costo) < distancia[grafo.G[actual.nodo][i].nodo])){
				distancia[grafo.G[actual.nodo][i].nodo] = actual.costo + grafo.G[actual.nodo][i].costo;
                previo[ grafo.G[actual.nodo][i].nodo ] = actual.nodo; 
				Q.push(Nodo(grafo.G[actual.nodo][i].nodo, actual.costo + grafo.G[actual.nodo][i].costo));
			}
		}
	}
	a = destino;
	camino( destino, previo);
	cout<<endl;
}
class TSP	
{
private:
	Grafo grafo;
	string iFile;
	string oFile;
	
	vector<int>odds;
	int **cost;
	
	vector<int> *adjList;
	void findOdds();
public:
	int n;
	int r;
	
	int **path_vals;
	
	int pathLength;
	
	vector<int> circuit;
	
	long int **graph;
  	vector<int>* adjlist;
	// Constructor
	TSP(string in, string out);
	// Destructor
	~TSP(){}
	void perfectMatching();
	void euler_tour(int start, vector<int> &path);
	void make_hamiltonian(vector<int> &path, int &pathCost);
	void findMST();
	int getMinIndex(int key[], bool mst[]);
	void make_shorter();
	void printResult();
	void printPath();
	int get_size(){return n;};
	void fillMatrix();
	int findBestPath(int start);
};
TSP::TSP(string in, string out){
	iFile = in;
	oFile = out;

	ifstream inStream;
	inStream.open("prueba.txt", ios::in);

	if(!inStream){
		cerr << "Can't open input file " << iFile << endl;
		exit(1);
	}
	int aux1, aux2,aux3;
	inStream >> aux1 >> aux2 >> aux3;
	
	n = aux1;
	grafo.V=n;
	grafo.A=aux2;
	r = aux3;
	graph = new long int*[n];
	for(int i = 0; i < n; i++){
		graph[i] = new long int[n];
		for(int j = 0; j < n; j++){
			graph[i][j] = 999;
		}
	}

	cost = new int*[n];
	for(int i = 0; i < n; i++){
		cost[i] = new int[n];
	}

	path_vals = new int*[n];
	for(int i = 0; i < n; i++){
		path_vals[i] = new int[n];
	}

	adjlist = new vector<int>[n];
    //READ DATA
	int c, x, y;
	while(!inStream.eof()){
		inStream >> c >> x >> y;
		graph[c-1][x-1]=graph[x-1][c-1]=y;
		grafo.G[c].push_back(Arista(x, y));
		grafo.G[x].push_back(Arista(c, y));
	}
	inStream.close();
}

void TSP::findMST() {

  int *key = new int[n];
  bool *included = new bool[n];
  int *parent = new int[n];

  for (int i = 0; i < n; i++) {
    key[i] = std::numeric_limits<int>::max();
    included[i] = false;
  }
  key[0] = 0;
  parent[0] = -1;

  for (int i = 0; i < n - 1; i++) {
    int k = getMinIndex(key, included);
    included[k] = true;
    for (int j = 0; j < n; j++) {
      if (graph[k][j] && included[j] == false && graph[k][j] < key[j]) {
          parent[j] = k;
          key[j] = graph[k][j];
      }
    }
  }
  for (int i = 0; i < n; i++) {
    int j = parent[i];
    if (j != -1) {
      adjlist[i].push_back(j);
      adjlist[j].push_back(i);
    }
  }
}

int TSP::getMinIndex(int key[], bool mst[]) {
  int min = std::numeric_limits<int>::max();
  int min_index;
  for (int i = 0; i < n; i++) {
    if (mst[i] == false && key[i] < min) {
      min = key[i];
      min_index = i;
    }
  }
  return min_index;
}



void TSP::findOdds() {

  for (int i = 0; i < n; i++) {
    if ((adjlist[i].size() % 2) != 0) {
      odds.push_back(i);
    }
  }
}


void TSP::perfectMatching() {
  int closest, length; //int d;
  std::vector<int>::iterator tmp, first;
  findOdds();
  while (!odds.empty()) {
    first = odds.begin();
    vector<int>::iterator it = odds.begin() + 1;
    vector<int>::iterator end = odds.end();
    length = std::numeric_limits<int>::max();
    for (; it != end; ++it) {
      
      if (graph[*first][*it] < length) {
        length = graph[*first][*it];
        closest = *it;
        tmp = it;
      }
    } 
    adjlist[*first].push_back(closest);
    adjlist[closest].push_back(*first);
    odds.erase(tmp);
    odds.erase(first);
  }
}

void TSP::euler_tour(int start, vector<int> &path){
	
	vector<int> *tempList = new vector<int>[n];
	for(int i = 0; i < n; i++){
		tempList[i].resize(adjlist[i].size());
		tempList[i] = adjlist[i];
	}

	stack<int> stack;
	int pos = start;
	path.push_back(start);
	while(!stack.empty() || tempList[pos].size() > 0){
		if(tempList[pos].empty()){
			path.push_back(pos);
			pos = stack.top();
			stack.pop();
		}
		else{
			stack.push(pos);
			int neighbor = tempList[pos].back();
			tempList[pos].pop_back();
			for(int i = 0; i < tempList[neighbor].size(); i++){
				if(tempList[neighbor][i] == pos){
					tempList[neighbor].erase(tempList[neighbor].begin()+i);
				}
			}
			pos = neighbor;
		}
	}
	path.push_back(pos);
}

void TSP::make_hamiltonian(vector<int> &path, int &pathCost){
	
	bool *visited = new bool[n];
	for(int i = 0; i < n; i++){
		visited[i] = 0;
	}
	
	pathCost = 0;

	int root = path.front();
	vector<int>::iterator cur = path.begin();
	vector<int>::iterator iter = path.begin()+1;
	visited[root] = 1;
	bool addMore = true;
	while(iter != path.end()){
		if(!visited[*iter]){
			pathCost += graph[*cur][*iter];
			cur = iter;
			visited[*cur] = 1;
			iter = cur + 1;
		}	
		else{
			iter = path.erase(iter);
		}
	}
	
	if ( iter != path.end() ){
		pathCost += graph[*cur][*iter];
	}
}

int TSP::findBestPath(int start){
	vector<int> path;
	euler_tour(start, path);
	int length;
	make_hamiltonian(path, length);
	return length;
}
void TSP::make_shorter(){
	twoOpt(graph, circuit, pathLength, n);
}

void TSP::printResult(){
  outputStream.open(oFile.c_str(), ios::out);
  outputStream << pathLength << endl;
  for (vector<int>::iterator it = circuit.begin(); it != circuit.end(); ++it) {
	if(graph[*it][*(it+1)]==999){
		cout<<*it+1<<" "<<*(it+1)+1<<" AQUI : ";
		dijkstraCAMINO(*it+1,*(it+1)+1,grafo);
	}
	else{
		outputStream << *it+1 << endl;
	}    
  }
  outputStream <<circuit.front()+1;
  outputStream.close();
};

void TSP::printPath(){
   int a =  0;
  for (vector<int>::iterator it = circuit.begin(); it != circuit.end()-1; ++it) {
	  if( graph[*it][*(it+1)] ==999){
		int aux = dijkstraCOSTO(*it+1,*(it+1)+1,grafo);
		a+= aux;
	  }
	  else{
		a+=graph[*it][*(it+1)] ;
	  }		
  }
  pathLength =a;
};
int main() {
	string input, output;
    output = "salida1";
	output.append(".txt");
	TSP tsp(input, output);
	cout << "tsp created" << endl;
	int tsp_size = tsp.get_size();
	tsp.findMST();
	cout << "MST created" << endl;
	tsp.perfectMatching();
	cout << "Matching completed" << endl;
	int best = INT_MAX;
	int bestIndex;
	for (long t = 0; t < tsp_size; t++) {
		int result = tsp.findBestPath(t);
		tsp.path_vals[t][0] = t; // set start
		tsp.path_vals[t][1] = result; // set end
		if (tsp.path_vals[t][1] < best) {
			bestIndex = tsp.path_vals[t][0];
			best = tsp.path_vals[t][1];
		}
	}
	tsp.euler_tour(bestIndex,tsp.circuit);
	tsp.make_hamiltonian(tsp.circuit,tsp.pathLength);
	tsp.make_shorter();
	tsp.make_shorter();
	tsp.make_shorter();
	tsp.make_shorter();
	tsp.make_shorter();
	cout << "Final length: " << tsp.pathLength << endl;
	tsp.printPath();
	tsp.printResult();
}