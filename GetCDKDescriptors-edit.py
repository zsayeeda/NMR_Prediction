
""" generated source for module GetCDKDescriptors """
from rdkit import Chem
import os
import numpy as np
from scipy.spatial import distance


# 
#  * ApplyCDKDescriptors.py
#  * Purpose: Calculate  descriptors from input SDF files
#  * Calling the constructor executes the algorithm
#  *
#  * 
class GetCDKDescriptors(object):
    """ generated source for class GetCDKDescriptors """
    ENGINE = DescriptorEngine(DescriptorEngine.ATOMIC)

    # 
    #   * Example main
    #   *
    #   
    @classmethod
    def main(cls, args):
        """ generated source for method main """
        inpath = "HMDB00001.sdf"
        outpath = "HMDB00001.csv"
        getDescriptorCSV(inpath, outpath, "")

    # 
    #   * Constructor, executing the algorithm
    #   *
    #   * @param: string The path to the input SDF file
    #   * @param: string The path to the output CSV file
    #   
    def __init__(self, inpath, outpath, descNamesStr):
        """ generated source for method __init__ """
        getDescriptorCSV(inpath, outpath, descNamesStr)

    # 
    #   * Calculate descriptors. Omits IPMolecularLearningDescriptor
    #   *
    #   * @param string path to SDF input file
    #   * @param string path to CSV output file
    #   * @param string comma-seperated list of descriptor names (if empty, all descriptors will be calculated)
    #   
    @classmethod
    def getDescriptorCSV(cls, sdfInputPath, csvOutputPath, descNamesStr):
        """ generated source for method getDescriptorCSV """
        mols = readMolecules(sdfInputPath)
        System.err.println("read " + len(mols) + " compounds")
        descriptors = cls.ENGINE.getDescriptorInstances()
        System.err.println("found " + len(descriptors) + " descriptors")
        descNames = Arrays.asList(descNamesStr.split(","))
        colNames = ArrayList()
        values = ArrayList()
        for desc in descriptors:
            if isinstance(desc, (IPAtomicLearningDescriptor, )):
                continue 
            tname = tnamebits[len(tnamebits)]
            if (0 > len(descNamesStr)) and (not descNames.contains(tname)):
                continue 
            while len(colNamesArr):
                colNamesArr[idx] = tname + "-" + colNamesArr[idx]
                idx += 1
            colNames.addAll(Arrays.asList(colNamesArr))
            for mol in mols:
                while i < atomCount:
                    atoms.add(mol.getAtom(i))
                    i += 1
                values.addAll(computeListsAtomic(mol, atoms, desc))
        ncol = len(values)
        nrow = mols.get(0).getAtomCount()
        fstream = FileWriter(csvOutputPath)
        out = BufferedWriter(fstream)
        out.write("SMILES,")
        c = 0
        while c < ncol:
            if c != 0:
                out.write(",")
            out.write(colNames.get(c))
            c += 1
        out.write("\n")
        r = 0
        while r < nrow:
            out.write(smi + ",")
            while c < ncol:
                if c != 0:
                    out.write(",")
                out.write("" + values.get(c)[r])
                c += 1
            out.write("\n")
            r += 1
        out.flush()

   
    


    ####### my code startted here #####

    def getAtomicDescriptor(self, sdf, descNamesStr):
        mols = sdf
        print("read {} compounds".format(len(mols)))
        descriptors = ['GetAtomicNum','GetBonds','GetDegree','GetExplicitValence','GetFormalCharge','GetHybridization','GetImplicitValence','GetIsAromatic','GetIsotope','GetNoImplicit','GetNumExplicitHs','GetNumImplicitHs','GetNumRadicalElectrons','GetTotalDegree','GetTotalNumHs','IsInRing']
        values = [[]]
        i = 0
            for mol in mols:
                atoms = mol.GetAtoms()
                for atom in atoms:
                    atomic_num = atom.GetAtomicNum()
                    bonds = len(atom.GetBonds())
                    get_degree = atom.GetDegree()
                    exp_valence = atom.GetExplicitValence()
                    formal_charge = atom.GetFormalCharge()
                    hyberdization = 0
                    hyber = atom.GetHybridization().name
                    if hyber == "SP":
                        hyberdization = 2
                    elif hyber == "SP2":
                        hyberdization = 3
                    elif hyber == "SP3":
                        hyberdization = 4
                    impli_valence = atom.GetImplicitValence()
                    bool_aromatic = atom.GetIsAromatic()
                    if bool_aromatic == "True":
                        is_aromatic = 1
                    else: 
                        is_aromatic = 0
                    get_iso = atom.GetIsotope()
                    impli = atom.GetNoImplicit()
                    num_of_implicit = 0
                    if impli == "True":
                        num_of_implicit = 1
                    else:
                        num_of_implicit = 0
                    num_of_exp_h = atom.GetNumExplicitHs()
                    num_of_impli_h = atom.GetNumImplicitHs()
                    num_of_radical_elec = atom.GetNumRadicalElectrons()
                    total_degree = atom.GetTotalDegree()
                    total_num_h = atom.GetTotalNumHs()
                    ring = atom.IsInRing()
                    is_in_ring = 0
                    if ring == "True":
                        is_in_ring = 1
                    else:
                        is_in_ring = 0
                    symbol = atom.GetSmarts().upper()
                    values.append([symbol,atomic_num,bonds,get_degree,exp_valence,formal_charge,hyberdization,impli_valence,is_aromatic,get_iso,num_of_implicit,num_of_exp_h,num_of_impli_h,num_of_radical_elec,total_degree,total_num_h,is_in_ring])
        return values.pop(0)

    #  Find all relevant hydrogen atoms in molecule
   
    def getHydrogenAtoms(self, sdf):
        hydrogen_positions = []
        mols = sdf
        if len(mols) != 0:
            for mol in mols:
                for atom in mol.GetAtoms():
                    if ((atom.GetSmarts()).upper() != "O") and ((atom.GetSmarts()).upper() != "N"):
                        hydrogen_positions.append(atom.GetIdx())
        return hydrogen_positions

   """ generated source for method getNearestAtoms """
    def getNearestAtoms(self, sdf):
        mols = sdf
        atom_distances = [[]]
        for mol in mols:
            atoms = mol.GetAtoms()
            i = 0
            for atom1 in atoms:
                distances = []
                j = 0
                for atom2 in atoms:
                    if i == j:
                        distances.append(99999.12)
                    else:
                        conf = mol.GetConformer()
                        atom1_point = conf.GetAtomPosition(i)
                        atom2_point = conf.GetAtomPosition(j)
                        firstPoint = (atom1_point.x,atom1_point.y,atom1_point.z)
                        secondPoint = (atom2_point.x,atom2_point.y,atom2_point.z)
                        dst = distance.euclidean(firstPoint, secondPoint)
                        distances.append(dst)
                    j=+1
                indices =[]
                d = distances[:]
                d.sort
                d_list = distances
                for x in range(len(distances)):
                    index = d_list.index(d[x])
                atom_distances.append(indices)
                i=+1
         return atom_distances.pop(0)


    @classmethod
    def getHoseCodesForMolecule(cls, mol):
        """ generated source for method getHoseCodesForMolecule """
        hoseG = HOSECodeGenerator()
        hoseCodes = ArrayList()
        atomCount = mol.getAtomCount()
        atoms = ArrayList()
        i = 0
        while i < atomCount:
            try:
                hoseCodes.add(hose)
                print "HOSE = " + hose + "\n"
            except CDKException as e:
                e.printStackTrace()
            i += 1
        return hoseCodes

    @classmethod
    def getSmiles(cls, m):
        """ generated source for method getSmiles """
        props = m.getProperties()
        for key in props.keySet():
            if key.__str__() == "STRUCTURE_SMILES" or key.__str__() == "SMILES":
                return props.get(key).__str__()
        g = SmilesGenerator()
        return g.createSMILES(m)

    @classmethod
    def computeLists(cls, mols, desc):
        """ generated source for method computeLists """
        print "computing descriptor " + getName(desc)
        values = computeDescriptors(mols, desc)
        return values

    @classmethod
    def computeListsAtomic(cls, mol, atoms, desc):
        """ generated source for method computeListsAtomic """
        print "computing descriptor " + getName(desc)
        values = computeDescriptorsAtomic(mol, atoms, desc)
        return values

    @classmethod
    def readMolecules(cls, filepath):
        """ generated source for method readMolecules """
        mols = Vector()
        file_ = File(filepath)
        if not file_.exists():
            raise IllegalArgumentException("file not found: " + filepath)
        list_ = List()
        try:
            if reader == None:
                raise IllegalArgumentException("Could not determine input file type")
            list_ = ChemFileManipulator.getAllAtomContainers(content)
            reader.close()
        except Exception as e:
            e.printStackTrace()
            return None
        for iAtomContainer in list_:
            mol = AtomContainerManipulator.removeHydrogens(mol)
            try:
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol)
            except Exception as e:
                e.printStackTrace()
            try:
                CDKHueckelAromaticityDetector.detectAromaticity(mol)
            except Exception as e:
                e.printStackTrace()
            if mol.getAtomCount() == 0:
                System.err.println("molecule has no atoms")
            else:
                mols.add(mol)
        return mols


    def readMoleculesString(self, sdf):
        """ generated source for method readMoleculesString """
        mols = Vector()
        if sdf == "":
            raise IllegalArgumentException("No sdf found" + sdf)
        list_ = List()
        try:
            if reader == None:
                raise IllegalArgumentException("Could not determine input file type")
            list_ = ChemFileManipulator.getAllAtomContainers(content)
            reader.close()
        except Exception as e:
            e.printStackTrace()
            return None
        for iAtomContainer in list_:
            try:
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol)
            except Exception as e:
                e.printStackTrace()
            try:
                CDKHueckelAromaticityDetector.detectAromaticity(mol)
            except Exception as e:
                e.printStackTrace()
            if mol.getAtomCount() == 0:
                System.err.println("molecule has no atoms")
            else:
                mols.add(mol)
        return mols

    @classmethod
    def computeDescriptors(cls, mols, descriptor):
        """ generated source for method computeDescriptors """
        vv = ArrayList()
        j = 0
        while j < getSize(descriptor):
            vv.add([None]*len(mols))
            j += 1
        i = 0
        while i < len(mols):
            if mols.get(i).getAtomCount() == 0:
                while j < getSize(descriptor):
                    vv.get(j)[i] = None
                    j += 1
            else:
                try:
                    if isinstance(res, (IntegerResult, )):
                        vv.get(0)[i] = float((res).intValue())
                    elif isinstance(res, (DoubleResult, )):
                        vv.get(0)[i] = (res).doubleValue()
                    elif isinstance(res, (DoubleArrayResult, )):
                        while j < getSize(descriptor):
                            vv.get(j)[i] = (res).get(j)
                            j += 1
                    elif isinstance(res, (IntegerArrayResult, )):
                        while j < getSize(descriptor):
                            vv.get(j)[i] = float((res).get(j))
                            j += 1
                    else:
                        raise IllegalStateException("Unknown idescriptor result value for '" + descriptor + "' : " + res.__class__)
                except Throwable as e:
                    System.err.println("Could not compute cdk feature " + descriptor)
                    e.printStackTrace()
                    while j < getSize(descriptor):
                        vv.get(j)[i] = None
                        j += 1
            while j < getSize(descriptor):
                if vv.get(j)[i] != None and (vv.get(j)[i].isNaN() or vv.get(j)[i].isInfinite()):
                    vv.get(j)[i] = None
                j += 1
            i += 1
        return vv

    @classmethod
    def computeDescriptorsAtomic(cls, mol, atoms, descriptor):
        """ generated source for method computeDescriptorsAtomic """
        vv = ArrayList()
        vv.add([None]*len(atoms))
        i = 0
        while i < len(atoms):
            if atoms.get(i) == None:
                vv.get(0)[i] = None
            else:
                try:
                    if isinstance(res, (IntegerResult, )):
                        vv.get(0)[i] = float((res).intValue())
                    elif isinstance(res, (DoubleResult, )):
                        vv.get(0)[i] = (res).doubleValue()
                    elif isinstance(res, (DoubleArrayResult, )):
                        vv.get(0)[i] = (res).get(0)
                    elif isinstance(res, (IntegerArrayResult, )):
                        vv.get(0)[i] = float((res).get(0))
                    else:
                        raise IllegalStateException("Unknown idescriptor result value for '" + descriptor + "' : " + res.__class__)
                except Throwable as e:
                    System.err.println("Could not compute cdk feature " + descriptor)
                    e.printStackTrace()
                    vv.get(0)[i] = 0.0
            if vv.get(0)[i] != None and (vv.get(0)[i].isNaN() or vv.get(0)[i].isInfinite()):
                vv.get(0)[i] = 0.0
            i += 1
        return vv

    @classmethod
    def getSize(cls, descriptor):
        """ generated source for method getSize """
        r = descriptor.getDescriptorResultType()
        if isinstance(r, (DoubleArrayResultType, )):
            return len(length)
        elif isinstance(r, (IntegerArrayResultType, )):
            return len(length)
        else:
            return 1

    @classmethod
    def getName(cls, descriptor):
        """ generated source for method getName """
        try:
            if name != None:
                return name
            else:
                return ""
        except Throwable as e:
            return ""


if __name__ == '__main__':
    import sys
    GetCDKDescriptors.main(sys.argv)
