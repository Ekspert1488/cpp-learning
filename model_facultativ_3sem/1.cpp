#include <iostream>
#include <fstream>

  

int main(){

 std:: cout << "Hello, World!" << std::endl;
 int n = 0;

std:: cout << " Print number of numbers:" << std::endl;
std::cin >> n;
std::ofstream out;

out.open("Hello_text.txt");

if (out.is_open()){

std:: cout << "Hello, World!" << std::endl;

for (int i = 1; i<=n; i++){

out <<' ' <<i;

 }
 }
    out.close();
    std::cout << "File has been created";
}