#include<memory>
#include<vector>
#include<algorithm>
class Line
{
public:
    Line();
    ~Line();
    Line(size_t size);

    Line(Line const & );
    Line(Line       &&);
    Line & operator=(Line const & );
    Line & operator=(Line       &&);
    
    size_t size() const;
    float & x(size_t it);
    float & y(size_t it);
    void bound_chk(size_t size);
private:
    std::vector<float> ptrx;
    std::vector<float> ptry;
    // Member data.
}; /* end class Line */

