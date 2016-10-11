asap_emt_driver.o: asap_emt_driver.cpp \
 Asap.h \
 AsapNamespace.h Debug.h asap_kim_api.h KimAtoms.h KimAsapPython.h Vec.h \
 IVec.h NeighborLocator.h AsapPython.h AsapObject.h Templates.h \
 KimTemplates.h Exception.h Atoms.h asap_emt_driver.h EMT.h SymTensor.h \
 Potential.h EMTParameterProvider.h TinyMatrix.h KimParameterProvider.h \
 EMTDefaultParameterProvider.h
asap_kim_api.o: asap_kim_api.cpp KimAsapPython.h \
 asap_kim_api.h KimAtoms.h Asap.h AsapNamespace.h Vec.h IVec.h \
 NeighborLocator.h AsapPython.h AsapObject.h Templates.h KimTemplates.h \
 Exception.h Atoms.h Potential.h SymTensor.h NeighborCellLocator.h \
 KimNeighborMIOPBCH.h KimNeighborLocator.h KimNeighborMIOPBCF.h \
 KimNeighborNEIGHPUREH.h KimNeighborNEIGHPUREF.h KimNeighborNEIGHRVECH.h \
 KimNeighborNEIGHRVECF.h Debug.h
AsapObject.o: AsapObject.cpp Asap.h AsapNamespace.h AsapObject.h
EMT.o: EMT.cpp EMT.h AsapPython.h KimAsapPython.h \
 Asap.h \
 AsapNamespace.h Atoms.h KimAtoms.h Vec.h IVec.h SymTensor.h Potential.h \
 Exception.h AsapObject.h EMTParameterProvider.h TinyMatrix.h \
 NeighborLocator.h Templates.h KimTemplates.h Timing.h mass.h Debug.h
EMTDefaultParameterProvider.o: EMTDefaultParameterProvider.cpp \
 EMTDefaultParameterProvider.h EMTParameterProvider.h AsapPython.h \
 KimAsapPython.h \
 Asap.h \
 AsapNamespace.h AsapObject.h TinyMatrix.h Exception.h
Exception.o: Exception.cpp AsapPython.h KimAsapPython.h \
 Exception.h Asap.h AsapNamespace.h Debug.h
KimAtoms.o: KimAtoms.cpp KimAtoms.h KimAsapPython.h \
 Asap.h \
 AsapNamespace.h Vec.h IVec.h Exception.h Debug.h
KimNeighborLocator.o: KimNeighborLocator.cpp KimNeighborLocator.h \
 NeighborLocator.h AsapPython.h KimAsapPython.h \
 AsapObject.h AsapNamespace.h Templates.h KimTemplates.h Exception.h \
 Asap.h Atoms.h KimAtoms.h Vec.h IVec.h Debug.h
KimNeighborMIOPBCF.o: KimNeighborMIOPBCF.cpp KimNeighborMIOPBCF.h \
 KimNeighborMIOPBCH.h KimNeighborLocator.h NeighborLocator.h AsapPython.h \
 KimAsapPython.h \
 AsapObject.h AsapNamespace.h Templates.h KimTemplates.h Exception.h \
 Asap.h Atoms.h KimAtoms.h Vec.h IVec.h
KimNeighborMIOPBCH.o: KimNeighborMIOPBCH.cpp KimNeighborMIOPBCH.h \
 KimNeighborLocator.h NeighborLocator.h AsapPython.h KimAsapPython.h \
 AsapObject.h AsapNamespace.h Templates.h KimTemplates.h Exception.h \
 Asap.h Atoms.h KimAtoms.h Vec.h IVec.h
KimNeighborNEIGHPUREF.o: KimNeighborNEIGHPUREF.cpp \
 KimNeighborNEIGHPUREF.h KimNeighborNEIGHPUREH.h KimNeighborLocator.h \
 NeighborLocator.h AsapPython.h KimAsapPython.h \
 AsapObject.h AsapNamespace.h Templates.h KimTemplates.h Exception.h \
 Asap.h Atoms.h KimAtoms.h Vec.h IVec.h
KimNeighborNEIGHPUREH.o: KimNeighborNEIGHPUREH.cpp \
 KimNeighborNEIGHPUREH.h KimNeighborLocator.h NeighborLocator.h \
 AsapPython.h KimAsapPython.h \
 AsapObject.h AsapNamespace.h Templates.h KimTemplates.h Exception.h \
 Asap.h Atoms.h KimAtoms.h Vec.h IVec.h
KimNeighborNEIGHRVECF.o: KimNeighborNEIGHRVECF.cpp \
 KimNeighborNEIGHRVECF.h KimNeighborNEIGHRVECH.h KimNeighborLocator.h \
 NeighborLocator.h AsapPython.h KimAsapPython.h \
 AsapObject.h AsapNamespace.h Templates.h KimTemplates.h Exception.h \
 Asap.h Atoms.h KimAtoms.h Vec.h IVec.h
KimNeighborNEIGHRVECH.o: KimNeighborNEIGHRVECH.cpp \
 KimNeighborNEIGHRVECH.h KimNeighborLocator.h NeighborLocator.h \
 AsapPython.h KimAsapPython.h \
 AsapObject.h AsapNamespace.h Templates.h KimTemplates.h Exception.h \
 Asap.h Atoms.h KimAtoms.h Vec.h IVec.h Debug.h
KimParameterProvider.o: KimParameterProvider.cpp KimParameterProvider.h \
 EMTDefaultParameterProvider.h EMTParameterProvider.h AsapPython.h \
 KimAsapPython.h \
 Asap.h \
 AsapNamespace.h AsapObject.h TinyMatrix.h Exception.h Debug.h
Matrix3x3.o: Matrix3x3.cpp Matrix3x3.h Vec.h Asap.h AsapNamespace.h
NeighborCellLocator.o: NeighborCellLocator.cpp NeighborCellLocator.h \
 AsapPython.h KimAsapPython.h \
 Asap.h \
 AsapNamespace.h Vec.h NeighborLocator.h AsapObject.h Templates.h \
 KimTemplates.h Exception.h Atoms.h KimAtoms.h IVec.h Matrix3x3.h \
 Timing.h Debug.h
NeighborLocatorInterface.o: NeighborLocatorInterface.cpp KimAsapPython.h \
 NeighborCellLocator.h AsapPython.h Asap.h AsapNamespace.h Vec.h \
 NeighborLocator.h AsapObject.h Templates.h KimTemplates.h Exception.h \
 Atoms.h KimAtoms.h IVec.h
Potential.o: Potential.cpp AsapPython.h KimAsapPython.h \
 Asap.h \
 AsapNamespace.h Atoms.h KimAtoms.h Vec.h IVec.h Potential.h Exception.h \
 AsapObject.h SymTensor.h Debug.h
Timing.o: Timing.cpp AsapPython.h KimAsapPython.h \
 Timing.h \
 Asap.h AsapNamespace.h TimingResults.h
Vec.o: Vec.cpp Vec.h Asap.h AsapNamespace.h
