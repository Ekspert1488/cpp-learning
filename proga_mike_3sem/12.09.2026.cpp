#include <iostream>
#include <string>
class Vehicle {
    protected:
        std:: string make ;
        std:: string model;
        int year;
        int price;
    public:
        Vehicle( std:: string make , std:: string model,int year,int price){
            this->make = make;
            this->model = model;
            this->price = price;
            this->year = year;
        }
        virtual void print(){
            std::cout << make <<"\t" <<  model <<"\t" <<  year <<"\t" <<  price <<  std ::endl;
         }



};
class Car: public virtual Vehicle {
    private:
        int num_doors;
        std::string body_style;
    public:
        Car(std:: string make , std:: string model,int year,int price , int num_doors , std::string body_style):
        Vehicle (make,model,year,price){
            this->num_doors = num_doors;
            this->body_style = body_style;
        }
        std::string operator[](unsigned index) const{
             switch (index)
        {
        case 0 : return make;
        case 1: return model;
        case 2: return std:: to_string (year) ;
        case 3: return  std:: to_string (price);
        case 4: return std::to_string(num_doors);
        case 5: return body_style;
        default: return "Bad Index";
            }   
        }
       

        



};
class Truck: public virtual  Vehicle{
    private:
         const char* bed_length;
         const char*towing_capacity;
    public:
        Truck(std:: string make , std:: string model ,int year,int price , const char* bed_length , const char* towing_capacity):
        Vehicle (make,model,year,price){
            this->bed_length = bed_length;
            this->towing_capacity = towing_capacity;

        }
        virtual void print(){
            std::cout << make <<" \t" <<  model <<"\t" <<  year <<"\t" <<  price << "\t" << bed_length << "\t" << towing_capacity << std ::endl;
        }

};







int main(){
    Car car{"Toyota", "Camry", 2022, 2900000, 4, "Sedan"};
    std::cout << car[0] << "\n" << car[1] << "\n" << car[2] << "\n" << car[3]<< std::endl;    
}   