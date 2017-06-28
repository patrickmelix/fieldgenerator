#!/usr/bin/python3
"""
See Documentation for explenations.
Patrick Melix
chemistry@melix.me
2017/04

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
__all__ = ['FieldGenerator']

import pickle
import sys, os
from ase import Atoms, io, neighborlist
import numpy as np
from functools import reduce
from operator import add,sub
import simpleeval


#########################
# Classes
#########################
class FieldGenerator():
   """ Class to store and create information about the complete ForceField and write it to a FIELD file """

   def __init__(self, systemName='systemName', workDir='', systemXYZFileName='xyz.reference',
                  xyzFileNames=['orig.xyz'], ffFileNames=['force_field.dat'], nMols=[1], molNames=['MoleculeName'], 
                  dataFileName='data.pickle',
                  unit='kcal', levcfg=0, imcon=2, engcfg=0.0,
                  floatFormat='{:20.6f}', fieldFileName='FIELD',
                  configFileName='CONFIG', disablePrinting=False):

      """ Create an instance and set the default values if not given by arguments """

      #general information on the entire system
      self.systemName = systemName
      self.workDir = workDir
      self.systemXYZFileName = systemXYZFileName
      #moolecular information
      if not len(set([len(xyzFileNames),len(ffFileNames), len(nMols)])) == 1:
         sys.exit("Length of lists for xyzFileNames, ffFileNames and nMols must be equal")
      self.xyzFileNames = xyzFileNames
      self.FFFileNames = ffFileNames
      self.nMols = nMols
      self.molNames = molNames
      #information for the config file, see DL_POLY manual section CONFIG-file
      self.levcfg = levcfg
      self.imcon = imcon
      self.engcfg = engcfg
      #processing information
      self.dataFileName = dataFileName
      self.unit = unit
      self.floatFormat = floatFormat
      self.fieldFileName = fieldFileName
      self.configFileName = configFileName
      self.disablePrinting = disablePrinting
      #status variables, do not change
      self._readInput = False
      self._builtStruc = False
      self._hasField = False
      #external functions
      self.conv = _SimpleEvalAllowNumbers() #allow for mathematical expressions inside the input file
   """END __init__"""

   def __getstate__(self):
      """ helper to enable pickling with self.conv which is not pickable """
      return {k: v for k, v in self.__dict__.items() if _shouldPickle(v)}

   def __setstate__(self, state):
      """ help recreate self.conv after loading from pickle """
      self.__dict__.update(state)
      self.conv = _SimpleEvalAllowNumbers() #allow for mathematical expressions inside the input file


   #########################
   # Public Functions to call
   #########################

   def create(self):
      """read all information and create all files"""
      self.readInput()
      self.buildStructure()

      #write CONFIG file
      self.writeConfig()

      #create data for field file
      self.makeField()

      #write the FIELD file
      self.writeField()
   """END create"""



   def updateLFSE(self, lfseList, iMol=0):
      """update the LFSE part of the field file
      lfseList should be a list of floats"""
      if iMol > len(self.nMols)-1:
         sys.exit("Molecule index out of range, only "+str(len(self.nMols))+" exist.")
      #delete the old part of fieldOut
      lines = len(self.fieldOut)
      i = 0
      while i < lines:
         i += 1
         if self.fieldOut[-1] == 'lfse':
            del self.fieldOut[-1]
            break
         del self.fieldOut[-1]
      #update with list given
      self._updateLFSEParamsFromList(lfseList,iMol)

      self._makeLFSESection()

      #write the FIELD file
      self.writeField()
   """END updateLFSE"""



   def updateFF(self, valueList, key, iMol=0):
      """update any part of the classical FF
      iMol is the number of the molecule we want to change (starts with 0), ordering as in creation
      valueList should be a list of floats in the desired ordering
      (same as in the input force_field.dat but without atom names)
      length of the list must be equal to the numbers of values in the inputFF"""
      if iMol > len(self.nMols)-1:
         sys.exit("Molecule index out of range, only "+str(len(self.nMols))+" exist.")
      #update the FF
      self._updateFFFromList(iMol, valueList,key)

      #update the fieldOut
      if key == 'ATOMS':
         self._makeAtomSection(iMol, update=True)
      elif key == 'BONDS':
         self._makeBondSection(iMol, update=True)
      elif key == 'ANGLES':
         self._makeAngleSection(iMol, update=True)
      elif key == 'DIHEDRALS':
         self._makeDihedralSection(iMol, update=True)
      elif key == 'IMPROPER':
         self._makeDihedralSection(iMol, update=True)
      elif key == 'INVERSION':
         self._makeInversionSection(iMol, update=True)
      elif key == 'VDW':
         self._makeVdWSection(update=True)
      else:
         sys.exit("Unknown key in updateFF")

      #write the new FIELD file
      self.writeField()
   """END updateFF"""





   def readInput(self):
      """read all input"""
      #read xyz
      self.atoms = []
      for i,filePath in enumerate(self.xyzFileNames):
         xyzFilePath = os.path.join(self.workDir,filePath)
         self.atoms.append(io.read(xyzFilePath))
         if not self.disablePrinting:
            print('Number of Atoms in Molecule '+str(i)+': ' + str(len(self.atoms[i])))

      #read atom names and all positions from xyzReferenceFile
      self._readAtomNames()
      
      #read FF
      self._readFF()
      self._checkInputFF()

      #set check to true
      self._readInput = True
   """END readInput"""




   def buildStructure(self):
      """build all lists defining the system"""
      if not self._readInput:
         sys.exit("You need to run self.readInput() first")
      #build the bondList
      self._makeBondList()
      if not self.disablePrinting:
         print('Number of Bonds: ' + ",".join([str(i) for i in self.nBonds]))

      #get angleList
      self._makeAngleList()
      if not self.disablePrinting:
         print('Number of Angles: ' + ",".join([ str(len(i)) for i in self.angleList ]))

      #get dihedralList
      self._makeDihedralList()
      if not self.disablePrinting:
         print('Number of Dihedral Angles: ' + ",".join([ str(len(i)) for i in self.dihedralList ]))

      #make optional sections
      self._makeImproperList()
      self._makeInversionList()

      #set status
      self._builtStruc = True
   """END buildStructure"""



   def save(self):
      """Saves instance to pickle file"""
      dataFilePath = os.path.join(self.workDir,self.dataFileName)
      tmpFileWrite = open(dataFilePath, 'wb')
      #dump to tmpFile
      pickle.dump(self, tmpFileWrite)
      tmpFileWrite.close()
   """END save"""



   @staticmethod
   def load(path):
      with open(path, 'rb') as f:
         return pickle.load(f)
   """END load"""






   def makeField(self):
      """produce content of FIELD file"""
      if not self._readInput:
         sys.exit("You need to run self.readInput() first")
      if not self._builtStruc:
         sys.exit("You need to run self.buildStructure() first")
      if not hasattr(self, 'fieldOut'):
         self.fieldOut = []
      else:
         del self.fieldOut[:]

      self.fieldOut.append(self.systemName)
      self.fieldOut.append('UNITS ' + self.unit)
      self.fieldOut.append('MOLECULES ' + str(len(self.nMols)))
      #intramolecular definitions
      for iMol in range(0,len(self.nMols)):
         self.fieldOut.append(self.molNames[iMol])
         self.fieldOut.append('NUMMOLS ' + str(self.nMols[iMol]))
         self._makeAtomSection(iMol)
         self._makeBondSection(iMol)
         self._makeAngleSection(iMol)
         self._makeDihedralSection(iMol)
         self._makeInversionSection(iMol)
         self.fieldOut.append('FINISH')
      #intermolecular definitions
      self._makeVdWSection()
      self.fieldOut.append('CLOSE')
      self._makeLFSESection()
      
      #set status
      self._hasField = True
   """END makeField"""




   def writeField(self):
      """write the FIELD file"""
      if not self._readInput:
         sys.exit("You need to run self.readInput() first")
      if not self._builtStruc:
         sys.exit("You need to run self.buildStructure() first")
      if not self._hasField:
         sys.exit("You need to run self.makeField() first")
      fieldFilePath = os.path.join(self.workDir,self.fieldFileName)
      write = open(fieldFilePath, 'w')
      for line in self.fieldOut:
         write.write(line)
         write.write('\n')
      write.close()
   """END writeField"""



   def writeConfig(self):
      """write the CONFIG file"""
      if not self._readInput:
         sys.exit("You need to run self.readInput() first")
      nAtoms = 0
      for i,mol in enumerate(self.atoms):
         nAtoms += len(mol)*self.nMols[i]
      configFilePath = os.path.join(self.workDir,self.configFileName)
      write = open(configFilePath, 'w')
      #systemname and box are static
      box = self.box
      write.write('{:>80}\n'.format(self.systemName))
      write.write((('{:10d}'*3)+'{:20f}\n').format(self.levcfg,self.imcon,nAtoms,self.engcfg))
      #if a periodic system is given (imcon > 0) write the box
      if self.imcon > 0:
         write.write((('{:20f}'*3)+'\n').format(*box[0]))
         write.write((('{:20f}'*3)+'\n').format(*box[1]))
         write.write((('{:20f}'*3)+'\n').format(*box[2]))

      #now all atoms
      xyzReferenceFilePath = os.path.join(self.workDir,self.systemXYZFileName)
      if not os.path.isfile(xyzReferenceFilePath):
         sys.exit('File ' + xyzReferenceFilePath + ' does not exist!')
      refFile = open(xyzReferenceFilePath, 'r')
      if self.imcon > 0:
         for i in range(0,3):
            next(refFile)
      i = 0
      for line in range(0,nAtoms):
         i += 1
         line = refFile.readline().split()
         #convert to float for proper alignement
         line[1:] = [ float(i) for i in line[1:] ]
         #allow for additional values for velocity and force
         nLines = int((len(line)-1)/3)
         write.write(('{:8}'+'{:10d}\n').format(line[0],i))
         for j in range(0,nLines):
            write.write(('{:20f}'*3+'\n').format(*line[(j*3)+1:((j+1)*3)+1]))
      write.close()
      refFile.close()
   """END writeConfig"""





   #########################
   # Private Functions
   #########################



   def _makeBondList(self):
      """Create a bondlist from the atoms
      if self.imcon > 0 we have a periodic structure"""
      periodic = ( self.imcon > 0 )
      #special case 2D slab
      twoD = (self.imcon == 6)
      if twoD:
         dimensions = 2
      elif periodic:
         dimensions = 3
      box = self.box
      #initiate bondList
      self.bondList = []
      self.uniqueBondList = []
      self.nBonds = []
      for i,molecule in enumerate(self.atoms):
         #make a copy of the original atoms
         tmpMol = molecule.copy()
         self.bondList.append({})
         nAtoms = len(tmpMol)
         if periodic:
            tmpMol.set_cell(self.box)
            if twoD:
               tmpMol.set_pbc((True,True,False))
            else:
               tmpMol.set_pbc((True,True,True))
         #create a list of cutoffs
         cutOff = []
         for j in range(0,nAtoms):
            cutOff.append(radii[tmpMol[j].symbol])
         #initiate neighborlist
         neighborList = neighborlist.NeighborList(cutOff,self_interaction=False,bothways=True)
         neighborList.update(tmpMol)

         for j in range(0,nAtoms):
            array, trash = neighborList.get_neighbors(j)
            #sort only to make results more ordered
            array = list(array)
            array.sort()
            #do not add duplicates from images
            self.bondList[i][str(j)] = list(set(array))

         #get number of bonds and make a list of unique bonds (no double counting)
         self.uniqueBondList.append([])
         nBonds = 0
         for iAtom,bonded in self.bondList[i].items():
            for jAtom in bonded:
               if ([int(iAtom),jAtom] in self.uniqueBondList[i]) or ([jAtom,int(iAtom)] in self.uniqueBondList[i]):
                  continue
               self.uniqueBondList[i].append([int(iAtom),jAtom])
               nBonds += 1
         self.nBonds.append(nBonds)
   """END makeBondList"""



   def _makeAngleList(self):
      """make a list of Angles"""
      self.angleList = []
      for i,molecule in enumerate(self.atoms):
         self.angleList.append([])
         for j,bonds in self.bondList[i].items():
            for k in bonds:
               for l in self.bondList[i][str(k)]:
                  if not int(j) == l:
                     self.angleList[i].append([int(j),k,l])
         #delete duplicates
         for angle in self.angleList[i]:
            del self.angleList[i][self.angleList[i].index(list(reversed(angle)))]
   """END makeAngleList"""



   def _makeDihedralList(self):
      """Make a list of all Dihedrals"""
      self.dihedralList = []
      for idx,molecule in enumerate(self.atoms):
         self.dihedralList.append([])
         for i,bonds in self.bondList[idx].items():
            for j in bonds:
               for k in self.bondList[idx][str(j)]:
                  if int(i) == k:
                     continue
                  for l in self.bondList[idx][str(k)]:
                     if (not l == k) and (not l == int(i)) and (not l == j):
                        self.dihedralList[idx].append([int(i),j,k,l])
         """delete duplicates"""
         for dihedral in self.dihedralList[idx]:
            if list(reversed(dihedral)) in self.dihedralList[idx]:
               del self.dihedralList[idx][self.dihedralList[idx].index(list(reversed(dihedral)))]
   """END makeDihedralList"""



   def _makeImproperList(self):
      """Make a list of all Improper"""
      self.improperList = []
      self.nImpropers = []
      for iMol,molecule in enumerate(self.atoms):
         improperData = self.inputFF[iMol]['IMPROPER']
         tmpData = {}
         self.improperList.append([])

         #create list with all atom numbers for each defined improper
         for i,improper in enumerate(improperData):
            tmpList = []
            for j,atomName in enumerate(improper[0:4]):
               tmpList.append([])
               #get all atoms with that name
               for iAtom,name in enumerate(self.names[iMol]):
                  if name == atomName:
                     tmpList[j].append(iAtom)
            #iterate over the four atom names
            for j,jList in enumerate(tmpList):
               idx = list(range(0,4))
               idx.remove(j)
               #iterate over atoms
               for jAtom in jList:
                  jBonds = self.bondList[iMol][str(jAtom)]
                  #for an improper angle only a three bonded atom can be the central atom
                  if not len(jBonds) == 3:
                     continue
                  tmpInv = []
                  #iterate over bonded atoms
                  for kAtom in jBonds:
                     for k in idx:
                        if kAtom in tmpList[k]:
                           tmpInv.append(kAtom)
                           break
                  if len(tmpInv) == 3:
                     tmpSet = [jAtom,*tmpInv]
                     exists = False
                     for prevSet in self.improperList[iMol]:
                        if set(tmpSet) == set(prevSet):
                           exists = True
                           break
                     if not exists:
                        self.improperList[iMol].append(tmpSet)
         self.nImpropers.append(len(self.improperList[iMol]))

      if not self.disablePrinting:
         print('Number of Improper Angles: ' + ",".join([str(i) for i in self.nImpropers]))
   """END makeImproperList"""



   def _makeInversionList(self):
      """Make a list of all Inverstion"""
      self.inversionList = []
      for iMol,molecule in enumerate(self.atoms):
         inversionData = self.inputFF[iMol]['INVERSION']
         tmpData = {}
         self.inversionList.append([])

         #create list with all atom numbers for each defined inversion
         for i,inversion in enumerate(inversionData):
            tmpList = []
            for j,atomName in enumerate(inversion[0:4]):
               tmpList.append([])
               #get all atoms with that name
               for iAtom,name in enumerate(self.names[iMol]):
                  if name == atomName:
                     tmpList[j].append(iAtom)
            #iterate over the four atom names
            for j,jList in enumerate(tmpList):
               idx = list(range(0,4))
               idx.remove(j)
               #iterate over atoms
               for jAtom in jList:
                  jBonds = self.bondList[iMol][str(jAtom)]
                  #for an inversion angle only a three bonded atom can be the central atom
                  if not len(jBonds) == 3:
                     continue
                  tmpInv = []
                  #iterate over bonded atoms
                  for kAtom in jBonds:
                     for k in idx:
                        if kAtom in tmpList[k]:
                           tmpInv.append(kAtom)
                           break
                  if len(tmpInv) == 3:
                     tmpSet = [jAtom,*tmpInv]
                     exists = False
                     for prevSet in self.inversionList[iMol]:
                        if set(tmpSet) == set(prevSet):
                           exists = True
                           break
                     if not exists:
                        self.inversionList[iMol].append(tmpSet)

      if not self.disablePrinting:
         print('Number of Inversion Angles: ' + ",".join([ str(len(i)) for i in self.inversionList ]))
   """END makeInversionList"""



   def _inCutOff(self,iMol,i,j):
      """determine if i and j are inside the cutoff distance"""
      dist = self.atoms[iMol].get_distance(i,j)
      iElement = self.atoms[iMol][i].symbol
      jElement = self.atoms[iMol][j].symbol
      cutOff = getCutOff(iElement,jElement)
      if cutOff > dist:
         return True
      else:
         return False
   """END inCutOff"""



   def _readAtomNames(self):
      """reads names for each atom and the simulation box"""
      self.names = []
      self.box = []
      xyzReferenceFilePath = os.path.join(self.workDir,self.systemXYZFileName)
      if not os.path.isfile(xyzReferenceFilePath):
         sys.exit('File ' + xyzReferenceFilePath + ' does not exist!')
      refFile = open(xyzReferenceFilePath, 'r')
      #box definition if we have a periodic system
      if self.imcon > 0:
         for i in range(0,3):
            line = refFile.readline()
            line = line.split()
            line = [ float(x) for x in line ]
            self.box.append(line)
         iLine = 3
      else:
         iLine = 0
      for iMol, nMol in enumerate(self.nMols): 
         self.names.append([])
         nLines = len(self.atoms[iMol])
         for i in range(0,nLines):
            lineList = refFile.readline().split()
            iLine += 1
            #check if atom name is already in use in another molecule if it appears the first time
            if not lineList[0] in self.names[iMol]:
               for j in range(0,iMol):
                  if lineList[0] in self.names[j]:
                     sys.exit("Atom names must be unique over all molecules")
            self.names[iMol].append(lineList[0])
         #fast forward to the next molecule
         for i in range(1,nMol):
            for j in range(0,nLines):
               iLine += 1
               next(refFile)
      #check that there is nothing more
      try:
         for line in refFile:
            if not line.strip() == '':
               sys.exit("More content in xyzReferenceFile then expected from the molecules given.")
      except:
         pass
      refFile.close()
   """END readAtomNames"""


   def _readFF(self):
      """reads ForceField from file"""
      self.inputFF = []
      for idx in range(0,len(self.nMols)):
         self.inputFF.append({})
         FFFilePath = os.path.join(self.workDir, self.FFFileNames[idx])
         if not os.path.isfile(FFFilePath):
            sys.exit('File ' + FFFilePath + ' does not exist!')
         readFile = open(FFFilePath, 'r')
         lines = []
         #read lines without comments or empty lines
         for line in readFile:
            if line[0] == '#':
               continue
            elif not line.strip():
               continue
            lines.append(line.split())
         i = -1
         iEnd = -1
         while i < len(lines)-1:
            i += 1
            #read key
            if i > iEnd:
               if not len(lines[i]) <= 3:
                  print(lines[i])
                  sys.exit('Error in input force-field! Section size is not right, check numbers of entries.')
               key = lines[i][0]
               if key == 'LFMM':
                  nItems = int(lines[i][1]) + int(lines[i][2])*6*2
               else:
                  nItems = int(lines[i][1])
               iEnd = i + nItems
               self.inputFF[idx][key] = []
               if key == 'LFMM':
                  self.inputFF[idx][key].append(lines[i][1:])
            else:
               self.inputFF[idx][key].append(lines[i])
         readFile.close()
   """END readFF"""


   def _updateFFFromList(self, iMol, valueList, key):
      """update a part of the force-field with the list given"""
      if key not in self.inputFF[iMol]:
         sys.exit('Unknown key given to update FF')

      #data to update
      data = self.inputFF[iMol][key]

      #have to have the same amount of values, potential and names not included!
      #if that is the case, update the content
      iStart = 0
      for i, line in enumerate(data):
         #find the first float
         for j, item in enumerate(line):
            if isFloat(item):
               break

         iEnd = iStart + len(line) - j
         line[j:] = valueList[iStart:iEnd]
         iStart = iEnd

      if not iStart == len(valueList):
         sys.exit('Length of list given does not match data in storage.')
   """END _updateFFFromList"""


   def _checkInputFF(self):
      """check the inputFF for errors like double entries"""
      for iMol in range(0,len(self.nMols)):
         atomData = self.inputFF[iMol]['ATOMS']
         #check that atomData only has unique entries, otherwise soething is wrong
         names = []
         for atom in atomData:
            names.append(atom[0])
         if not len(set(names)) == len(names):
            sys.exit("Non-unique label for atom in inputFF")
   """END checkInputFF"""





   def _makeAtomSection(self, iMol, update=False):
      """make the atom section of for the FIELD file"""
      section = []
      nAtoms = len(self.names[iMol])
      section.append('ATOMS ' + str(nAtoms))
      atomData = self.inputFF[iMol]['ATOMS']
      tmpData = {}
      for i in range(0,nAtoms):
         name = self.names[iMol][i]
         if name not in tmpData:
            tmpData[name] = findInFF([name],atomData)
         mass = float(self.conv.eval(tmpData[name][1]))
         charge = float(self.conv.eval(tmpData[name][2]))
         section.append(name + '      ' + (self.floatFormat*2).format(mass,charge))
      #add to or update the field list
      if update:
         overwrite(self.fieldOut,iMol,section)
      else:
         self.fieldOut.extend(section)
   """END makeAtomSection"""


   def _makeBondSection(self, iMol, update=False):
      section = []
      nBonds = self.nBonds[iMol]
      section.append('BONDS ' + str(nBonds))
      bondList = self.bondList[iMol]
      bondData = self.inputFF[iMol]['BONDS']
      tmpData = {}

      for bond in self.uniqueBondList[iMol]:
         bondName = [ self.names[iMol][n] for n in bond ]
         bondName.sort()
         bondNameKey = ' '.join(bondName)
         if bondNameKey not in tmpData:
            tmpData[bondNameKey] = findInFF(bondName,bondData)
         potential = tmpData[bondNameKey][2]
         values = [float(self.conv.eval(k)) for k in tmpData[bondNameKey][3:]]
         tmpBond = [ i+1 for i in bond ]
         section.append((('{:8} '*3)+(self.floatFormat*len(values))).format(potential,*tmpBond,*values))
      #add to or update the field list
      if update:
         overwrite(self.fieldOut,iMol,section)
      else:
         self.fieldOut.extend(section)
   """END makeBondSection"""



   def _makeAngleSection(self, iMol, update=False):
      """make the angle section of for the FIELD file"""
      section = []
      nAngles = len(self.angleList[iMol])
      section.append('ANGLES ' + str(nAngles))
      angleList = self.angleList[iMol]
      angleData = self.inputFF[iMol]['ANGLES']
      tmpData = {}

      for angle in angleList:
         angleName = [ self.names[iMol][n] for n in angle ]
         revAngleName = list(reversed(angleName))
         angleNameKey = ' '.join(angleName)
         revAngleNameKey = ' '.join(revAngleName)
         if angleNameKey not in tmpData:
            if revAngleNameKey not in tmpData:
               tmpData[angleNameKey] = findInFF(angleName,angleData)
            else:
               angleNameKey = revAngleNameKey
               angleName = revAngleName
         potential = tmpData[angleNameKey][3]
         values = [float(self.conv.eval(k)) for k in tmpData[angleNameKey][4:]]
         tmpAngle = [ i+1 for i in angle ]
         section.append((('{:8} '*4)+(self.floatFormat*len(values))).format(potential,*tmpAngle,*values))
      #add to or update the field list
      if update:
         overwrite(self.fieldOut,iMol,section)
      else:
         self.fieldOut.extend(section)
   """END makeAngleSection"""




   def _makeDihedralSection(self, iMol, update=False):
      """make the dihedral section of for the FIELD file"""
      section = []
      nDihedrals = len(self.dihedralList[iMol])
      nImproper = self.nImpropers[iMol]
      section.append('DIHEDRALS ' + str(nDihedrals+nImproper))
      dihedralList = self.dihedralList[iMol]
      dihedralData = self.inputFF[iMol]['DIHEDRALS']
      tmpData = {}

      for dihedral in dihedralList:
         dihedralName = [ self.names[iMol][n] for n in dihedral]
         revDihedralName = list(reversed(dihedralName))
         dihedralNameKey = ' '.join(dihedralName)
         revDihedralNameKey = ' '.join(revDihedralName)
         if dihedralNameKey not in tmpData:
            if revDihedralNameKey not in tmpData:
               tmpData[dihedralNameKey] = findInFF(dihedralName,dihedralData)
            else:
               dihedralNameKey = revDihedralNameKey
               dihedralName = revDihedralName
         potential = tmpData[dihedralNameKey][4]
         values = [float(self.conv.eval(k)) for k in tmpData[dihedralNameKey][5:]]
         tmpDihedral = [ i+1 for i in dihedral ]
         section.append((('{:8} '*5)+(self.floatFormat*len(values))).format(potential,*tmpDihedral,*values))

      #add the improper part
      self._makeImproperSection(iMol, section)
      #add to or update the field list
      if update:
         overwrite(self.fieldOut,iMol,section)
      else:
         self.fieldOut.extend(section)
   """END makeDihedralSection"""







   def _makeInversionSection(self, iMol, update=False):
      """make the Inversion Section"""
      if not 'INVERSION' in self.inputFF[iMol] or len(self.inputFF[iMol]['INVERSION']) == 0:
         self.fieldOut.append('INVERSION  0')
         return

      section = []
      inversionData = self.inputFF[iMol]['INVERSION']
      tmpData = {}

      #if the inversionList was not built before
      if not hasattr(self, 'inversionList'):
         self._makeInversionList()

      #print the inversions
      section.append('INVERSION ' + str(len(self.inversionList[iMol])))
      #iterate over each inversion
      for inversion in self.inversionList[iMol]:
         inversionName = [ self.names[iMol][n] for n in inversion]
         inversionNameKey = ' '.join(inversionName)
         #find set, since order is not relevant
         tmpData[inversionNameKey] = findSetInFF(inversionName,inversionData)
         potential = tmpData[inversionNameKey][4]
         values = [float(self.conv.eval(k)) for k in tmpData[inversionNameKey][5:]]
         tmpInversion = [ x+1 for x in inversion ]
         section.append((('{:8} '*5)+(self.floatFormat*len(values))).format(potential,*tmpInversion,*values))
      #add to or update the field list
      if update:
         overwrite(self.fieldOut,iMol,section)
      else:
         self.fieldOut.extend(section)
   """END makeInversionSection"""




   def _makeImproperSection(self, iMol, section):
      """make the Improper Section
      note that impropers are added to dihedrals, never call this function directly"""
      if not 'IMPROPER' in self.inputFF[iMol] or len(self.inputFF[iMol]['IMPROPER']) == 0:
         return

      improperData = self.inputFF[iMol]['IMPROPER']
      tmpData = {}
      
      #iterate over each improper
      for improper in self.improperList[iMol]:
         improperName = [ self.names[iMol][n] for n in improper]
         improperNameKey = ' '.join(improperName)
         #find set, since order is not relevant
         tmpData[improperNameKey] = findSetInFF(improperName,improperData)
         potential = tmpData[improperNameKey][4]
         values = [float(self.conv.eval(k)) for k in tmpData[improperNameKey][5:]]
         tmpImproper = [ x+1 for x in improper ]
         section.append((('{:8} '*5)+(self.floatFormat*len(values))).format(potential,*tmpImproper,*values))
   """END makeImproperSection"""




   def _makeCrosstermSection(self, iMol, update=False):
      """make the Crossterm Section"""
      pass
   """END makeCrosstermSection"""





   def _makeVdWSection(self, update=False):
      """make the VdW Section"""
      #gather information from different molecules
      nAtoms = 0
      VdWData = []
      for iMol, mol in enumerate(self.atoms):
         nAtoms += len(self.inputFF[iMol]['ATOMS'])
         VdWData.extend(self.inputFF[iMol]['VDW'])

      section = []
      nVdW = nAtoms
      for i in range(1,nAtoms):
         nVdW += i
      section.append('VDW ' + str(nVdW))
      for i in range(0,nAtoms):
         iAtomName = VdWData[i][0]
         iAtomValues = VdWData[i][1:]
         for j in range(i,nAtoms):
            jAtomName = VdWData[j][0]
            jAtomValues = VdWData[j][1:]
            if not jAtomValues[0] == iAtomValues[0]:
               sys.exit('VdW-Potential missmatch')
            epsilon = float(self.conv.eval(iAtomValues[1])) + float(self.conv.eval(jAtomValues[1]))
            sigma = np.sqrt(float(self.conv.eval(iAtomValues[2])) * float(self.conv.eval(jAtomValues[2])))
            #lennart-jones 12-6
            if iAtomValues[0] == 'lj':
               values = [ sigma, epsilon ]
            #buckingham,  1.12246208 is 2^(1/6).  184000 and 2.25 constants from the MM3 form
            elif iAtomValues[0] == 'buck':
               values = [ 184000*epsilon, 2*sigma*1.22462048/12, 2.25*epsilon*(pow(2*sigma*1.122462048,6)) ]
            else:
               sys.exit('Requested VdW-Potential ' + iAtomValues[0] + ' not known.')

            section.append((('{:8} '*3)+(self.floatFormat*len(values))).format(iAtomName,jAtomName,iAtomValues[0],*values))
      #add to or update the field list
      if update:
         #only one VDW section exists, so we can just pass 0
         overwrite(self.fieldOut,0,section)
      else:
         self.fieldOut.extend(section)
   """END makeVdWSection"""






   def _makeLFSESection(self):
      """make the LFSE Section"""
      if not any(list( ('LFMM' in i) for i in self.inputFF)):
         return
      first = True
      for iMol in range(0,len(self.nMols)):
         #if no lfse in this molecule skip
         if (not 'LFMM' in self.inputFF[iMol]) or (len(self.inputFF[iMol]['LFMM']) == 0):
            continue
         if first:
            self.fieldOut.append('lfse')
            first = False

         lfseData = self.inputFF[iMol]['LFMM']
         nAtomTypes = int(lfseData[0][0])
         nBondTypes = int(lfseData[0][1])
         elements = []
         spinLabel = ['lo', 'hi']
         #read metal information
         for i in range(0,nAtomTypes):
            elements.append({})
            elements[i]['name'] = lfseData[i+1][0]
            elements[i]['nBonds'] = int(lfseData[i+1][1])
            elements[i]['lowspin'] = ' '.join(lfseData[i+1][2:7])
            elements[i]['highspin'] = ' '.join(lfseData[i+1][7:])
            if elements[i]['highspin'] == elements[i]['lowspin']:
               elements[i]['nSpinStates'] = 1
            else:
               elements[i]['nSpinStates'] = 2

         #gather info about all lfse centers
         for element in elements:
            atoms = []
            nAtoms = 0
            spinStates = element['nSpinStates']
            self.fieldOut.append('mc')
            metalName = element['name']
            for i,n in enumerate(self.names[iMol]):
               if n == metalName:
                  atom = _lfseAtom()
                  atom.index = i
                  atom.name = metalName
                  atom.nBonds = element['nBonds']
                  atom.lowspin = element['lowspin']
                  atom.highspin = element['highspin']
                  atom.bonds = self.bondList[iMol][str(i)]
                  if not len(atom.bonds) == atom.nBonds:
                     print(atom.bonds)
                     print(atom.nBonds)
                     sys.exit('Problem with number of bonds on metal-center ' + atom.name + ' with index ' + str(i+1))
                  atoms.append(atom)
                  nAtoms += 1

            #write bonding information
            for i in range(0,nAtoms):
               if i == 0:
                  string = 'atoms = (('
               else:
                  string = '         ( '
               bonding = [ atoms[i].index, *atoms[i].bonds ]
               bondedNames = [ self.names[iMol][x] for x in atoms[i].bonds ]
               #bonding order might not be the same for all centers due to periodic boundaries
               #therefore we must ensure that all bonding is in the same order with respect to the elements
               #reference will simply be the first bonding info
               if i == 0:
                  refBondingNames = bondedNames[:]
               #if not first item make now sure that the bonded atoms are in the same order with respect to their names
               elif bondedNames == refBondingNames:
                  pass
               #else we have to reorder them, so that the name pattern matches
               else:
                  oldBondList = atoms[i].bonds[:]
                  #create a list of tuples to help the sorting
                  tupleList = [ (oldBondList,bondedNames[i]) for i in range(0,len(oldBondList)) ]
                  #now sort tupleList so that the bondnames match refBondingNames
                  newBondList = []
                  for refName in refBondingNames:
                     #find the first occurance of refName in bondedNames
                     index = bondedNames.index(refName)
                     #add the corresponding number to the new list
                     newBondList.append(oldBondList[index])
                     #delete the occurance from the original lists
                     del oldBondList[index]
                     del bondedNames[index]
                  #make the new bonding list
                  bonding = [ atoms[i].index, *newBondList ]

               bonding = [ str(i+1) for i in bonding ]
               string += ', '.join(bonding)
               if i == nAtoms-1:
                  string += '))'
               else:
                  string += '),'
               self.fieldOut.append(string)

            for iSpin in range(0,spinStates):
               #write coefficients and save morse for later writing
               tmpData = {}
               #not really orbitals, but I had no better name... Call it contribution dimension if zou want.
               orbitals = []
               #spaces for alignement, no other reason
               orbitalLabels = ['as_'+spinLabel[iSpin]+' ', 'apx_'+spinLabel[iSpin], 'apy_'+spinLabel[iSpin], 'ads_'+spinLabel[iSpin], 'aee   ']
               morse = []
               #three values for the morse potential
               morse.append([])
               morse.append([])
               morse.append([])
               #number of contibutions is 5 for the low spin and 4 for the high spin (aee only once)
               nContrib = 5
               start = 0
               if iSpin == 1:
                  start = 5
                  nContrib = 4
               for iContrib in range(0,nContrib):
                  orbitals.append([])
                  metalName = element['name']
                  #all atoms of this type should have the same bonding
                  atomIdx = [ i for i,x in enumerate(atoms) if x.name == metalName ][0]
                  #for every bond print coeffs
                  iBond = 0
                  for idxBonded in atoms[atomIdx].bonds:
                     iBond += 1
                     bondedName = self.names[iMol][idxBonded]
                     bond = [metalName, bondedName]
                     bondKey = ''.join(bond)
                     if bondKey not in tmpData:
                        tmpData[bondKey] = findInFF(bond,lfseData,lfse=True)
                        if not len(tmpData[bondKey]) == spinStates*6 - (spinStates-1):
                           sys.exit('Number of entries in the LFSE part of the force-field does not match the required number')
                     #only one morse-potential per bond
                     if iContrib == 0:
                        back = -1
                        if spinStates == 2 and iSpin == 0:
                           back = -2
                        morse[0].append(float(self.conv.eval(tmpData[bondKey][back][2])))
                        morse[1].append(float(self.conv.eval(tmpData[bondKey][back][3])))
                        morse[2].append(float(self.conv.eval(tmpData[bondKey][back][4])))
                     values = tmpData[bondKey][iContrib+start][2:]
                     values = [ float(self.conv.eval(f)) for f in values ]
                     if iBond == 1:
                        string = orbitalLabels[iContrib] + ' = (('
                     else:
                        string = '          ('
                     string += (((self.floatFormat+', ')*6)+self.floatFormat).format(*values)
                     if iBond == len(atoms[atomIdx].bonds):
                        string += ' ))'
                     else:
                        string += ' ),'
                     #append this line
                     self.fieldOut.append(string)

                  #empty line between orbitals
                  self.fieldOut.append('')

               #write morse
               morseTitles = ['e0', 'r0', 'k']
               for i,item in enumerate(morse):
                  string = 'morse_' + morseTitles[i] + '_' + spinLabel[iSpin] + ' = ( '
                  string += (((self.floatFormat+', ')*(len(item)-1))+self.floatFormat).format(*item)
                  string += ' )'
                  self.fieldOut.append(string)
               self.fieldOut.append('')
            #end spin-loop

            #print spin configurations
            self.fieldOut.append('lo_spin = ' + element['lowspin'])
            self.fieldOut.append('hi_spin = ' + element['highspin'])

            #end mc section
            self.fieldOut.append('end mc')
            self.fieldOut.append('')
         #end molecule loop
      #end lfse section
      self.fieldOut.append('')
      self.fieldOut.append('! default is off')
      self.fieldOut.append('! flip_spins = off')
      self.fieldOut.append('! flip_spins = on')
      self.fieldOut.append('end lfse')
   """END makeLFSESection"""




   def _updateLFSEParamsFromList(self, inList, iMol):
      """read new LFSEParameters from list"""
      lfmmData = self.inputFF[iMol]['LFMM']
      nMetals = int(lfmmData[0][0])
      itemInListStart = 0
      for i in range(nMetals+1,len(lfmmData)):
         overwrite = lfmmData[i]
         itemInListEnd = itemInListStart + len(overwrite[2:])
         overwrite[2:] = inList[itemInListStart:itemInListEnd]
         itemInListStart = itemInListEnd
   """END updateLFSEParamsFromList"""


#########################
# General Functions
#########################

def findInFF(searchList, data, lfse=False):
   """find a search list in the provided list of lists"""
   items = len(searchList)
   lfseList = []
   for dataList in data:
      if (dataList[0:items] == searchList) or (dataList[0:items] == list(reversed(searchList))):
         if not lfse:
            return dataList
         else:
            lfseList.append(dataList)
   if lfse and len(lfseList) > 0:
      return lfseList
   sys.exit('Item ' + str(searchList) + ' not found in Force Field!')
"""END findInFF"""



def findSetInFF(searchList, data):
   """find a search list in the provided list of lists"""
   items = len(searchList)
   for dataList in data:
      if set(dataList[0:items]) == set(searchList):
         return dataList
   sys.exit('Set ' + str(searchList) + ' not found in Force Field!')
"""END findSetInFF"""


def getCutOff(iElement,jElement):
   """get the CutOff for different Element-Combinations"""
   global radii
   global radiiFactor
   try:
      return (radii[iElement]+radii[jElement])*radiiFactor
   except:
      sys.exit("Could not get cutoff for elements "+iElement+" "+jElement)
"""END getCutOff"""


def isFloat(item):
   try:
      float(item)
      return True
   except ValueError:
      return False



def overwrite(oldList,iMol,newSection):
   key = newSection[0]
   #search for the first line to replace
   found = False
   iFound = 0
   for i, item in enumerate(oldList):
      if item == key:
         #we are looking for the iMol-th Molecule
         if iFound == iMol:
            found = True
            break
         else:
            iFound += 1
   #if line is not found
   if not found:
      sys.exit("Key not found in old list "+str(key)+", cannot overwrite old data!")
   #set start and end positions in oldList
   iStart = i
   iEnd = iStart + len(newSection)
   #compare lengths of sections
   if len(oldList[iEnd].split()) > 3:
      sys.exit("Length of new data does not match old data")

   oldList[iStart:iEnd] = newSection[:]
"""END overwrite"""


def _shouldPickle(v):
   """ Check if this should be pickled """
   if isinstance(v, _SimpleEvalAllowNumbers):
      return False
   else:
      return True

#########################
# Global constants
#########################
global radii
global radiiFactor
radiiFactor = 1.1
radii = {}
radii[ 'H'] = 0.30 
radii['He'] = 0.99 
radii['Li'] = 1.52 
radii['Be'] = 1.12 
radii[ 'B'] = 0.88 
radii[ 'C'] = 0.77 
radii[ 'N'] = 0.70 
radii[ 'O'] = 0.66 
radii[ 'F'] = 0.64 
radii['Ne'] = 1.60 
radii['Na'] = 1.86 
radii['Mg'] = 1.60 
radii['Al'] = 1.43 
radii['Si'] = 1.17 
radii[ 'P'] = 1.10 
radii[ 'S'] = 1.04 
radii['Cl'] = 0.99 
radii['Ar'] = 1.92 
radii[ 'K'] = 2.31 
radii['Ca'] = 1.97 
radii['Sc'] = 1.60 
radii['Ti'] = 1.46 
radii[ 'V'] = 1.31 
radii['Cr'] = 1.25 
radii['Mn'] = 1.29 
radii['Fe'] = 1.26 
radii['Co'] = 1.25 
radii['Ni'] = 1.24 
radii['Cu'] = 1.28 
radii['Zn'] = 1.33 
radii['Ga'] = 1.41 
radii['Ge'] = 1.22 
radii['As'] = 1.21 
radii['Se'] = 1.17 
radii['Br'] = 1.14 
radii['Kr'] = 1.97 
radii['Rb'] = 2.44 
radii['Sr'] = 2.15 
radii[ 'Y'] = 1.80 
radii['Zr'] = 1.57 
radii['Nb'] = 1.41 
radii['Mo'] = 1.36 
radii['Tc'] = 1.35 
radii['Ru'] = 1.33 
radii['Rh'] = 1.34 
radii['Pd'] = 1.38 
radii['Ag'] = 1.44 
radii['Cd'] = 1.49 
radii['In'] = 1.66 
radii['Sn'] = 1.62 
radii['Sb'] = 1.41 
radii['Te'] = 1.37 
radii[ 'I'] = 1.33 
radii['Xe'] = 2.17 
radii['Cs'] = 2.62 
radii['Ba'] = 2.17 
radii['La'] = 1.88 
radii['Ce'] = 1.818
radii['Pr'] = 1.824
radii['Nd'] = 1.814
radii['Pm'] = 1.834
radii['Sm'] = 1.804
radii['Eu'] = 2.084
radii['Gd'] = 1.804
radii['Tb'] = 1.773
radii['Dy'] = 1.781
radii['Ho'] = 1.762
radii['Er'] = 1.761
radii['Tm'] = 1.759
radii['Yb'] = 1.922
radii['Lu'] = 1.738
radii['Hf'] = 1.57 
radii['Ta'] = 1.43 
radii[ 'W'] = 1.37 
radii['Re'] = 1.37 
radii['Os'] = 1.34 
radii['Ir'] = 1.35 
radii['Pt'] = 1.38 
radii['Au'] = 1.44 
radii['Hg'] = 1.52 
radii['Tl'] = 1.71 
radii['Pb'] = 1.75 
radii['Bi'] = 1.70 
radii['Po'] = 1.40 
radii['At'] = 1.40 
radii['Rn'] = 2.40 
radii['Fr'] = 2.70 
radii['Ra'] = 2.20 
radii['Ac'] = 2.00 
radii['Th'] = 1.79 
radii['Pa'] = 1.63 
radii[ 'U'] = 1.56 
radii['Np'] = 1.55 
radii['Pu'] = 1.59 
radii['Am'] = 1.73 
radii['Cm'] = 1.74 
radii['Bk'] = 1.70 
radii['Cf'] = 1.86 
radii['Es'] = 1.86 
radii['Fm'] = 2.00 
radii['Md'] = 2.00 
radii['No'] = 2.00 
radii['Lr'] = 2.00 
radii['Rf'] = 2.00 
radii['Db'] = 2.00 
radii['Sg'] = 2.00 
radii['Bh'] = 2.00 
radii['Hs'] = 2.00 
radii['Mt'] = 2.00 
radii['Ds'] = 2.00 
radii['Rg'] = 2.00 
radii['Cn'] = 2.00 
radii['Nh'] = 2.00 
radii['Fl'] = 2.00 
radii['Mc'] = 2.00 
radii['Lv'] = 2.00 
radii['Ts'] = 2.00 
radii['Og'] = 2.00 


#########################
# Private Classes
#########################
class _SimpleEvalAllowNumbers(simpleeval.SimpleEval):
   """Subclass SimpleEval to allow numbers to be passed to it"""
   def eval(self, expr):
      if isinstance(expr,float):
         tmp = str(expr)
      elif isinstance(expr,int):
         tmp = str(expr)
      else:
         tmp = expr
      return super(_SimpleEvalAllowNumbers, self).eval(tmp)



class _lfseAtom():
   def __init__(self):
      self.index = -1
      self.name = ''
      self.nBonds = 0
      self.bonds = []
      self.lowspin = ''
      self.highspin = ''
   """END __init__"""
"""END lfseAtom"""

if __name__ == "__main__":
   instance = FieldCreator()
   instance.create()
