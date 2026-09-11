#include <iostream>


class Node {
        
    public:
        int val;
        Node* next;
        Node ( int value, Node* nextNode = nullptr):
            val(value),next(nextNode)
        {}

        


};



void addFive (Node* node) {
    node->val+=5;
}
void print(Node* n1){
    std::cout << n1->val << n1->next->val << n1->next->next-> val << n1->next->next->next-> val  << std::endl;

}
void print_val(const Node& node){
    std:: cout << node.val << std:: endl;
}

int main() {
   Node* n1 = new Node(1);
   Node* n2 = new Node(2);
   Node* n3 = new Node(3);
   Node* n4 = new Node(4);
   n1->val = 10;
   n2->val = 20;
   n3->val = 30;
   n4->val = 40;
   n1->next = n2;
   n2->next =n3;
   n3->next = n4;
   print(n1);
   addFive(n1->next->next);
   print(n1);
   print_val(*n1);
   delete n4;
   delete n3;
   delete n2;
   delete n1;











}