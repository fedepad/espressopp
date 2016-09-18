#include <boost/python.hpp>
#include <vector>
#include <set>
#include <string>

      template <typename T>
        std::vector<T> python_list_to_vector (
                    const boost::python::object& obj
                    )
        {
                std::vector<T> vect(boost::python::len(obj));
                    for (unsigned long i = 0; i < vect.size(); ++i)
                            {
                                        vect[i] = boost::python::extract<T>(obj[i]);
                            }
                        return vect;
        }


    // if values in a list are unique we can save the content in a set which allows for logN search (std::vector instead should be sorted...NO!)
    template <typename T>
        std::set<T> python_list_to_set (
                    const boost::python::object& obj
                    )
        {
                //std::set<T> setn(boost::python::len(obj));
                std::set<T> setn;
                    //for (unsigned long i = 0; i < setn.size(); ++i)
                    for (unsigned long i = 0; i < boost::python::len(obj); ++i)
                            {
                                        //vect[i] = boost::python::extract<T>(obj[i]);
                                        setn.insert(boost::python::extract<T>(obj[i]));
                            }
                        return setn;
        }


   template <typename T>
       boost::python::list vector_to_python_list (
                   const std::vector<T>& vect
                   )

       {
               boost::python::list obj;
                   for (unsigned long i = 0; i < vect.size(); ++i)
                               obj.append(vect[i]);
                       return obj;
       }
