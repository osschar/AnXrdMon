#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedenum;
#pragma link C++ nestedtypedef;

#pragma link C++ class SRange+;

#pragma link C++ class SXrdFileInfo+;
#pragma link C++ class SXrdUserInfo+;
#pragma link C++ class SXrdServerInfo+;

#pragma link C++ class    SXrdReq+;
#pragma link C++ class    std::vector<SXrdReq>;
#pragma link C++ typedef vSXrdReq_t;

#pragma link C++ class    SXrdIoInfo+;

#pragma link C++ function dump_xrdfar_tree_ioinfo;

#endif
