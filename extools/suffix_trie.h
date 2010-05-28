/*
 * suffix_trie.h
 *
 *  Created on: May 17, 2010
 *      Author: Vlad
 */

#ifndef SUFFIX_TRIE_H_
#define SUFFIX_TRIE_H_

#include <string>
#include <map>
#include <vector>

template <class T>
class Node {
	public:
		std::map<T, Node> edge_list_;
		Node(){};
		Node(T w)
			: label_ (w){};
		  int InsertEdge(std::vector<T> p, int start, int end);
		  bool CheckEdge(std::vector<T> p);

		static int countNode_;
		T label_;
		int id_;

};

bool DEBUG = false;

template <class T>
int Node<T>::InsertEdge(std::vector<T> p, int start, int end){
  if (DEBUG) std::cout << "s " << start << " " << end << "\n||Insert||\n";
/*
 * Always start at root node at work forward through FSA
 */
	Node* currNode = this;
	typename std::map<T, Node>::iterator it;
	std::map<T, Node> *curr_edge_list_;
	std::pair<typename std::map<T ,Node>::iterator,bool> checkEdge;


/*
 * Start at the begining of the start and go to the last word
 */
	for(int i=start;i<= end; i++ )
	{
		curr_edge_list_ = &currNode->edge_list_;

		if (DEBUG) std::cout << "Word: " << p[i] << "\n";
		Node n_(p[i]);
		//check to see if node is already defined
	

		checkEdge = curr_edge_list_->insert(std::pair<T,Node> (p[i], n_));
		if(checkEdge.second==false)
		{
		  if(DEBUG)
		    {		
		      std::cout << "element " << currNode->label_ << " already has "  << checkEdge.first->second.label_ << std::endl;
		    }
		}
		//if(!curr_edge_list_->count(p[i]))

		if (DEBUG) {
		  std::cout << currNode->label_ << " map contains:\n";	       
		  for ( it=curr_edge_list_->begin() ; it != curr_edge_list_->end(); it++ )
		    std::cout << (*it).first << " => " << (*it).second.label_ << std::endl;
		  std::cout <<"\n";
		}
		currNode = &checkEdge.first->second;

	}
	if(DEBUG){std::cout <<"\n";}
return 1;
}


/*
 * Check to see if a path exists
 */
template <class T>
bool Node<T>::CheckEdge(std::vector<T> p){

	/*
 * Always start at root node at work forward through FSA
 */
  if (DEBUG)std::cout << "check:";
  Node* currNode = this;
  typename std::map<T, Node>::iterator it;
  std::map<T, Node> *curr_edge_list_;

  bool exists = true;
  int i =0;
/*
 * Start at the begining of the start and go to the last word
 */
  while(exists)
    {
      curr_edge_list_ = &currNode->edge_list_;
      
      if (DEBUG)      std::cout << "Word: " << p[i] << "\n";
      it = curr_edge_list_->find(p[i]);
      
      if(it != curr_edge_list_->end() )
	{
	  if (DEBUG)std::cout << "exists";
	  currNode= &it->second;
	}
      else
	{
	  if (DEBUG)  std::cout<< "false";
	  exists = false;
	  
	}
      
      i++;
    }
  
  if (p.size() == i-1)
    return true;
  else
    return false;
  //std::cout << std::endl << p.size() <<" Matched " << i-1 << "\n";

}


#endif /* SUFFIX_TRIE_H_ */
