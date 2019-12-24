#include <vector>
#include <iostream>
#include "Line.h"	


Line::Line(){};
Line::Line(Line const& other) {
	for (size_t i=0; i<other.data.size();i++){
		data.push_back(other.data[i]);
	}
};
Line::Line(Line&& other){
	for (size_t i=0; i<other.data.size();i++){
		data.push_back(other.data[i]);
	}
	other.data.clear();
};
Line& Line::operator=(Line const& other){
	if(this == &other){
		return *this;
	}else{
		this->data = other.data;
		return *this;
	}
};
Line& Line::operator=(Line&& other){
	if(this == &other){
		return *this;
	}else{
		this->data = other.data;
		other.data.clear();
		return *this;
	}
};
Line::Line(size_t size) {
	for(size_t i=0; i<size; i++){
		data.push_back(std::pair<float, float> (0, 0));
	}
// # TODO: maybe use reserve is the better way to to this thing,
//         but the code below is not working.
//	data.reserve(size+1);
//	std::cout << "size: " << size  << std::endl;
//	std::cout << "capacity: " << data.capacity() << std::endl;
};
Line::~Line(){};
size_t Line::size() const{
	return this->data.size();
};
float const & Line::x(size_t it) const{
	return this->data.at(it).first;
};
float& Line::x(size_t it){
	return this->data.at(it).first;
};
float const & Line::y(size_t it) const{
	return this->data.at(it).second;
};
float& Line::y(size_t it){
	return this->data.at(it).second;
};

