//
//  GenericFactory.hpp
//  FiniteElementAnalysis
//
//  Created by Thomas Kummer on 05.05.17.
//  Copyright © 2017 Thomas Kummer. All rights reserved.
//

#ifndef GenericFactory_hpp
#define GenericFactory_hpp

#include <stdio.h>
#include <exception>
#include <typeinfo>


namespace LifeV {

template<class factory_type>
class GenericFactory;

template<class factory_type>
class iCreator
{
public:
    virtual factory_type* Create(const std::string& key) const = 0;
    virtual ~iCreator() {}
};


template<class factory_type, class product_type>
class Creator : public iCreator<factory_type>
{
public:
    Creator(const std::string& key);
    virtual factory_type* Create(const std::string& key) const
    {
        return new product_type;
    }

};
    
template<class factory_type>
class GenericFactory
{
public:
    
    static GenericFactory& Instance()
    {
        static GenericFactory factory;
        return factory;
    }
    
    factory_type* Create(const std::string& key) const
    {
        auto i = _makers.find(key);
        if (i == _makers.end())
        {
            std::cerr << "Unrecognized object type!";
        }
        iCreator<factory_type>* maker = i->second;
        return maker->Create(key);
    }
    
    void RegisterMaker(const std::string& key, iCreator<factory_type>* maker)
    {
        if (_makers.find(key) != _makers.end())
        {
            std::cerr << "Multiple makers for given key!";
        }
        _makers[key] = maker;
    }
    
    void printRegisteredCreators(const std::string& type)
    {
        std::cout << "Registered " << type << " derived classes: " << std::endl;
        for (auto& creator : _makers) std::cout << " - " << creator.first << std::endl;
    }
    
private:
    
    std::map<std::string, iCreator<factory_type>*> _makers;

};

    
template<class factory_type, class product_type>
Creator<factory_type, product_type>::Creator(const std::string& key)
{
    GenericFactory<factory_type>::Instance().RegisterMaker(key, this);
}

    
#define REGISTER(factory_type, product_type) static Creator<factory_type, product_type> (creator_ ## product_type)(#product_type)
#define CREATE(factory_type, product_type) GenericFactory<factory_type>::Instance().Create(product_type)
#define PRINT_FACTORY(factory_type) GenericFactory<factory_type>::Instance().printRegisteredCreators(#factory_type)

}



#endif /* GenericFactory_hpp */
