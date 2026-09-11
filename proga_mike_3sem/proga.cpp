#include <iostream>
class Table 
{   
    friend void setPrice(Table&, int price);
    friend int getPrice (Table&);
    public: 
        Table (float length, float width, float height , int price);
        void print();
        void SetHeight(float p_height);
        float GetHeight();
        Table (const Table &t);
        Table &up (float height){   
            this->height += height;
            return *this;}
        void print_count ();
        ~Table (){
            --count;
            std::cout << " The table has been deleted " << std::endl;
        }
        
        
    private:
        
       
        static inline  unsigned int count{};
    protected:
        float width;
        float height;
        int price;
        float length;
       
     
};
Table::Table(float length, float width, float height, int price){
            count++;
            this->length = length;
            this->width = width;
            this->height = height;
            this->price = price;
            std:: cout << "The Table has been created" << std::endl ;
        }

void Table::print() 
{
    std::cout << "Length: " << this->length << " Width: " << this->width << " Height:" << this->height << std::endl;

}
void Table::SetHeight( float p_height){
    height = p_height;
}
float Table::GetHeight () {
    return height;
}
Table::Table (const Table &t){
    length = t.length;
    height = t.height + 0.1;
    width = t.width;
    ++count;

}
void setPrice(Table &table , int price){
    table.price = price;
}
int getPrice (Table &table){
    return table.price;
}
void Table::print_count(){
    std:: cout << "Number of created tables: " << count << std::endl;
}
class Monitor 
{
    protected:
        int number_of_devices;
    public:
        Monitor(int number_of_devices){
            this->number_of_devices = number_of_devices;
        }
        void print(){
            std::cout << "Number of devices: " << number_of_devices << std::endl;
        }
};



class DeskTable : public  virtual Table , public virtual Monitor
{
 private:
    float square;
 public:
 
    DeskTable (float length, float width, float height, int price , int number_of_devices): Table(length,width,height,price), Monitor(number_of_devices)
    {
        this->square = length * width;

    }
    void print_square(){
        std:: cout << "Square of this table:" << square << std::endl;
    }
    void print_length()
        {
            std:: cout << "Length of this table: " << length << std::endl;
        }
    
    void print (DeskTable &table ){
        table.Monitor::print();
        table.Table::print();
        

    }
    
    




};


int main()
{
    DeskTable table{1.5, 2.0, 0.75, 3000, 1};
    table.DeskTable::print(table);

   

}

