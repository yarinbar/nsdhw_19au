#include <iostream>
#include<vector>
#include<utility>
#include<stddef.h>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <iostream>
using namespace std;


class Line
{
public:
    Line(){
        this->m_size = 0;
    };
    Line(Line const &ref_line){
        this->m_size = ref_line.m_size;
        this->m_x = ref_line.m_x;
        this->m_y = ref_line.m_y;
    };
    /*Line(Line &&ref_line){
        *(this->m_x.begin()) = *(ref_line.m_x.begin());
        *(this->m_y.begin()) = *(ref_line.m_y.begin());
        (this->m_size) = (ref_line.m_size);
    }*/
    Line &operator=(Line const &ref_line){
        *(this->m_x.begin()) = *(ref_line.m_x.begin());
        *(this->m_y.begin()) = *(ref_line.m_y.begin());
        (this->m_size) = (ref_line.m_size);
    };
    //Line &operator=(Line &&);
    Line(size_t size){
        this->m_size = size;
        this->m_x.resize(size);
        this->m_y.resize(size);
    };
    ~Line() {};
    size_t size() const;
    
    
    size_t size(){
        return m_size;
    }
    float x(size_t it) const{
        return m_x[it];
    }
    
    float &x(size_t it){
        return *(m_x.data()+it);
    }    
    float y(size_t it) const{
        return m_y[it];
    }
    
    float &y(size_t it){
        return *(m_y.data()+it);
    }
    
private:
    // Member data.
    void check_range(size_t it) const;
    vector<float> m_x;
    vector<float> m_y;
    size_t m_size;
}; /* end class Line */ 

int main(int, char **)
{
    
    Line line(3);
    
    line.x(0) = 1; line.y(0) = 1;
    line.x(1) = 1; line.y(1) = 3;
    line.x(2) = 2; line.y(2) = 5;
    

    Line line2(line);
    line2.x(0) = 9;

    std::cout << "line: number of points = " << line.size() << std::endl;
    for (size_t it=0; it<line.size(); ++it)
    {
        std::cout << "point " << it << ":"
                  << " x = " << line.x(it)
                  << " y = " << line.y(it) << std::endl;
    }

    std::cout << "line2: number of points = " << line.size() << std::endl;
    for (size_t it=0; it<line.size(); ++it)
    {
        std::cout << "point " << it << ":"
                  << " x = " << line2.x(it)
                  << " y = " << line2.y(it) << std::endl;
    }
        
    return 0;
}