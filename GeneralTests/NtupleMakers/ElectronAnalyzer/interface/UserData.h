#ifndef CMS1_BaseTypes_h
#define CMS1_BaseTypes_h
// Description:     user data base class
// 
// Usage: 
// - event data - T
// - candidate data - vector<T> (vector of candidates each containing T).
// 
// Examples: 
// - vector of tracks in an event: T = vector<LorentzVector>
// - vector of tracks inside a candidate: T = vector<LorentzVector>, but the data is vector<T>
// - MET inside a candidate: T = float, but the data is vector<float>
//
// Original Author: Dmytro Kovalskyi
//
// $Author: dmytro $
// $Date: 2007/05/22 07:12:39 $
// $Revision: 1.2 $
//
#include "DataFormats/Math/interface/LorentzVector.h"
#include <string>
#include <memory>
#include <vector>

namespace cms1 {
   typedef math::XYZTLorentzVector LorentzVector;
   template <class T> class UserData
     {
      public:
	// If alias prefix is not set, no alias is created
	UserData(const std::string& name, 
		 const std::string& name_prefix, 
		 const std::string& alias_prefix, 
		 bool candidate = false)
	  {
	     isCandidate_ = candidate;
	     theName = name_prefix+name;
	     thePtrToData = &theData;
	     thePtrToDataVector = &theDataVector;
	     theAlias = alias_prefix;
	     if (theAlias.length()>0) theAlias += name;
	  }
	
	T*                   get()                      { return thePtrToData; }     // get access to data
	T**                  getAddress()               { return &thePtrToData; }     // ROOT needs something like that
	std::vector<T>*      getVector()                { return thePtrToDataVector; }     // get access to data
	std::vector<T>**     getVectorAddress()         { return &thePtrToDataVector; }     // ROOT needs something like that
	const std::string&   name()                     { return theName; }
	const std::string&   alias()                    { return theAlias; }
	bool                 isCandidate()              { return isCandidate_; }
	void                 addData( const T& data )   { theData = data; theDataVector.push_back( data ); } 
	void                 clearData()                { theData = T(); theDataVector.clear(); }
	
      private:
	// use of T and vector<T> allows for conditional unwrapping of data 
	// this is used in changing from event based ntuples to candidate based
	T theData;
	T* thePtrToData;
	std::vector<T> theDataVector;
	std::vector<T>* thePtrToDataVector;
	std::string theName;
	std::string theAlias;
	bool isCandidate_;
     };
   
   // Define some useful user data types
   typedef UserData<int>                                       UserDataInt;
   typedef UserData<std::vector<int> >                         UserDataInt1D;
   typedef UserData<float>                                     UserDataFloat;
   typedef UserData<std::vector<float> >                       UserDataFloat1D;
   typedef UserData<LorentzVector>                             UserDataP4;
   typedef UserData<std::vector<LorentzVector> >               UserDataP41D;
}
#endif
