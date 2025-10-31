//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.potential.groovy

import ffx.potential.MolecularAssembly
import ffx.potential.bonded.AminoAcidUtils
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Bond
import ffx.potential.bonded.MSNode
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.cli.PotentialScript
import ffx.potential.cli.SaveOptions
import ffx.potential.parameters.BioType
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static org.apache.commons.io.FileUtils.copyFile
import static org.apache.commons.io.FileUtils.createParentDirectories
import static org.apache.commons.io.FilenameUtils.*

/**
 * The SaveAsPDBNew script saves a file as a PDB file.
 * It mirrors SaveAsPDB behavior but is provided as a separate entry point.
 *
 * Usage:
 * ffxc SaveAsPDBNew [options] <filename>
 */
@Command(description = " Save the system as a PDB file.", name = "SaveAsPDBNew")
class SaveAsPDBNew extends PotentialScript {

  @Mixin
  SaveOptions saveOptions

  /**
   * --fs or ---firstSnapshot Provide the number of the first snapshot to be written.
   */
  @Option(names = ['--fs', '--firstSnapshot'], paramLabel = "-1", defaultValue = "-1",
      description = 'First snapshot to write out (indexed from 0).')
  private int firstSnapshot

  /**
   * --ls or ---lastSnapshot Provide the number of the last snapshot to be written.
   */
  @Option(names = ['--ls', '--lastSnapshot'], paramLabel = "-1", defaultValue = "-1",
      description = 'Last snapshot to write out (indexed from 0).')
  private int lastSnapshot

  /**
   * --si or --snapshotIncrement Provide the number of the snapshot increment.
   */
  @Option(names = ['--si', '--snapshotIncrement'], paramLabel = "1", defaultValue = "1",
      description = 'Increment between written snapshots.')
  private int snapshotIncrement

  /**
   * --wd or --writeToDirectories Provide the number of the snapshot increment.
   */
  @Option(names = ['--wd', '--writeToDirectories'], paramLabel = "false", defaultValue = "false",
      description = 'Write snapshots to numbered subdirectories.')
  private boolean writeToDirectories = false

  /**
   * --wp or --writeProperties Copy the property file to each subdirectory.
   */
  @Option(names = ['--cp', '--copyProperties'], paramLabel = "true", defaultValue = "true",
      description = 'Copy the property file to numbered subdirectories (ignored if not writing to subdirectories).')
  private boolean copyProperties = true

  /**
   * The final argument is an XYZ or ARC coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in XYZ or ARC format.')
  private String filename = null

  private Map<String, List<Integer>> moleculeAtomTypeDict = new HashMap<>()

  private Map<Integer, Atom> copiedAtomMap = new HashMap<>()

  /**
   * SaveAsPDBNew Constructor.
   */
  SaveAsPDBNew() {
    this(new groovy.lang.Binding())
  }

  /**
   * SaveAsPDBNew Constructor.
   * @param binding Groovy Binding to use.
   */
  SaveAsPDBNew(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  SaveAsPDBNew run() {

    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }
    // todo testing
    List<Atom> testAtoms = activeAssembly.getAtomList(true)
    for (Atom a : testAtoms) {
      logger.info(a.toString())
    }
    return this

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath()

    logger.info("\n Saving PDB for " + filename)

    // Use the current base directory, or update if necessary based on the given filename.
    String dirString = getBaseDirString(filename)

    String name = removeExtension(getName(filename)) + ".pdb"
    File saveFile = new File(dirString + name)

    SystemFilter openFilter = potentialFunctions.getFilter()

    // Handle snapshot writing to separate files.
    if (firstSnapshot >= 0) {
      openFilter.readNext(true)
      boolean resetPosition = true
      int counter = 0
      int snapshotCounter = 0
      logger.info(" Writing snapshots from " + firstSnapshot + " to " + lastSnapshot + " with increment " + snapshotIncrement)

      while (openFilter.readNext(resetPosition)) {
        resetPosition = false
        int offset = counter - firstSnapshot
        if (counter >= firstSnapshot && counter <= lastSnapshot && offset % snapshotIncrement == 0) {
          File snapshotFile
          if (writeToDirectories) {
            String subdirectory = concat(dirString, snapshotCounter.toString())
            snapshotFile = new File(concat(subdirectory, name))
            createParentDirectories(snapshotFile)
            if (copyProperties) {
              String propertyFile = activeAssembly.getProperties().getString("propertyFile")
              if (propertyFile != null) {
                String propertyFilename = getName(propertyFile)
                File copyOfPropFile = new File(concat(subdirectory, propertyFilename))
                logger.info("\n Copying properties to " + copyOfPropFile.toString())
                copyFile(new File(propertyFile), copyOfPropFile)
              }
            }
          } else {
            snapshotFile = new File(concat(dirString, removeExtension(name) + "." + counter.toString() + ".pdb"))
          }
          potentialFunctions.versionFile(snapshotFile)
          saveOptions.preSaveOperations(activeAssembly)
          MolecularAssembly newAssembly = buildNewAssembly()
          PDBFilter snapshotFilter = new PDBFilter(snapshotFile, newAssembly, newAssembly.getForceField(), newAssembly.getProperties())
          snapshotFile.append(activeAssembly.getCrystal().toCRYST1())
          snapshotFilter.writeFile(snapshotFile, true)
          snapshotCounter++
        }
        counter++
      }
      return this
    }

    saveFile = potentialFunctions.versionFile(saveFile)

    int numModels = openFilter.countNumModels()

    if (numModels == 1) {
      // Write one model and return.
      saveOptions.preSaveOperations(activeAssembly)
      MolecularAssembly newAssembly = buildNewAssembly()
      PDBFilter pdbFilter = new PDBFilter(saveFile, newAssembly, newAssembly.getForceField(), newAssembly.getProperties())
      saveFile.append(activeAssembly.getCrystal().toCRYST1())
      pdbFilter.writeFile(saveFile, true)
      return this
    }

    // Write out multiple models into a single PDB file.
    saveFile.write("MODEL        1\n")
    saveFile.append(activeAssembly.getCrystal().toCRYST1())
    MolecularAssembly firstAssembly = buildNewAssembly()
    saveOptions.preSaveOperations(activeAssembly)
    PDBFilter saveFilter = new PDBFilter(saveFile, firstAssembly, null, null)
    saveFilter.writeFile(saveFile, true, false, false)
    saveFile.append("ENDMDL\n")
    int modelNum = 1
    saveFilter.setModelNumbering(modelNum)


    // Iterate through the rest of the models in an arc or pdb.
    if (openFilter != null && (openFilter instanceof XYZFilter || openFilter instanceof PDBFilter)) {
      try {
        while (openFilter.readNext(false)) {
          saveOptions.preSaveOperations(activeAssembly)
          MolecularAssembly nextAssembly = buildNewAssembly()
//          PDBFilter modelFilter = new PDBFilter(saveFile, nextAssembly, nextAssembly.getForceField(), nextAssembly.getProperties())
          PDBFilter modelFilter = new PDBFilter(saveFile, nextAssembly, null, null)
          modelFilter.setModelNumbering(modelNum)
          modelFilter.writeFile(saveFile, true, true, false)
          modelNum++
        }
      } catch (Exception e) {
        // Do nothing.
      }
      // Add a final "END" record.
      saveFile.append("END\n")
    }
    return this
  }

  /**
   * Build a new MolecularAssembly from the current activeAssembly by identifying DNA chains,
   * waters, and ions, and constructing appropriate PDB-ready hierarchy.
   */
  private MolecularAssembly buildNewAssembly() {
    // Refresh the molecule atom type dictionary from the current force field.
    moleculeAtomTypeDict.clear()
    Map<String, BioType> bioTypes = activeAssembly.getForceField().getBioTypeMap()
    Map<String, ArrayList<BioType>> moleculeDict = new HashMap<>()
    if (!bioTypes.values().isEmpty()) {
      for (BioType bioType : bioTypes.values()) {
        if (!moleculeDict.containsKey(bioType.moleculeName)) {
          moleculeDict.put(bioType.moleculeName, new ArrayList<>())
          moleculeAtomTypeDict.put(bioType.moleculeName, new ArrayList<>())
        }
        ArrayList<BioType> moleculeBioTypes = moleculeDict.get(bioType.moleculeName)
        moleculeBioTypes.add(bioType)
        ArrayList<Integer> moleculeAtomTypes = moleculeAtomTypeDict.get(bioType.moleculeName)
        moleculeAtomTypes.add(bioType.atomType)
      }
    }

    List<Atom> atoms = new ArrayList<>(activeAssembly.getAtomList()) // make a copy

    // create a copied atom map
    copiedAtomMap.clear()
    for (Atom a : atoms) {
      Atom aCopy = new Atom(a.xyzIndex, a, a.getXYZ(new double[3]), -1, 'A' as char, "A")
      copiedAtomMap.put(a.xyzIndex as Integer, aCopy)
    }

    // adding bonds for the copied atoms
    for (Atom a : atoms) {
      int i = a.xyzIndex
      List<Atom> bondedAtoms = new ArrayList<>(a.get12List())
      Atom copy = copiedAtomMap.get(i as Integer)
      List<Atom> alreadyAddedAtoms = new ArrayList<>(copy.get12List())
      for (Atom ba : bondedAtoms) {
        int bi = ba.xyzIndex
        Atom bondedCopy = copiedAtomMap.get(bi as Integer)
        if (!alreadyAddedAtoms.contains(bondedCopy)) {
          Bond b = new Bond(copy, bondedCopy)
        }
      }
    }

    // todo - delete
//    atom    231     41  N     "N-Terminal NH3+"            7   14.0070  4
//    atom    232     42  H     "N-Terminal H3N+"            1    1.0080  1
//    atom    233     30  C     "C-Terminal COO-"            6   12.0110  3
//    atom    234     31  O     "C-Terminal COO-"            8   15.9990  1

    int fivePrimeHType = moleculeAtomTypeDict.get("5-Hydroxyl DNA").get(1) // 0 = O5* ; 1 = H5T
    int nTermNH3Type = moleculeAtomTypeDict.get("N-Terminal GLY").get(0) // 231 - its in all N-terminal amino acids, except proline
    int waterType = moleculeAtomTypeDict.get("Water").get(0) // 0 = O ; 1 = H
    int naType = moleculeAtomTypeDict.get("Sodium Ion").get(0)
    int clType = moleculeAtomTypeDict.get("Chloride Ion").get(0) // todo - could make list of all ion types
    List<Atom> fivePrimeHs = new ArrayList<>()
    List<Atom> nTermNH3s = new ArrayList<>()
    List<Atom> waterOs = new ArrayList<>()
    List<Atom> ions = new ArrayList<>()

    // get water oxygens, ions, and 5' oxygens of DNA
    for (Atom atom : atoms) {
      int at = atom.getType()
      if (at == waterType) { // get water O
        waterOs.add(atom)
      } else if (at == naType || at == clType) {
        ions.add(atom)
      } else if (at == fivePrimeHType) { // get HO5' (start of the DNA chain)
        fivePrimeHs.add(atom)
      } else if (at == nTermNH3Type) { // get N of NH3+ (start of protein chain)
        nTermNH3s.add(atom)
      }
    }

    String chainCharOptions = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" // maximum of 26 chains

    // Make all DNA chains
    int c = 0
    MSNode polymers = new MSNode()
    for (Atom atom : fivePrimeHs) {
      Character chain = chainCharOptions.charAt(c)
      Polymer polymer = addDNAChain(atom, chain)
      polymers.add(polymer)
      c++
    }

    // Make all protein chains
    for (Atom atom : nTermNH3s) {
      Character chain = chainCharOptions.charAt(c)
      Polymer polymer = addProteinChain(atom, chain)
      polymers.add(polymer)
      c++
    }

    // create new molecular assembly with polymers created above
    MolecularAssembly newAssembly = new MolecularAssembly("NewAssembly", polymers, activeAssembly.getForceField().getProperties())
    newAssembly.setForceField(activeAssembly.getForceField())

    // add the water molecules created above
    Character chain = 'A'
    int resNum = 1
    for (Atom atom : waterOs) {
      List<Atom> hAtoms = atom.get12List()
      if (hAtoms.size() != 2) {
        logger.severe("WATER OXYGEN DOES HAS MORE THAN TWO BONDS: " + atom.toString())
      }
      Atom h1 = hAtoms.get(0)
      h1.setName("H1")
      Atom h2 = hAtoms.get(1)
      h2.setName("H2")
      Atom[] water = [atom, h1, h2]
      for (Atom a : water) {
        a.setHetero(true)
        a.setResName("HOH")
        a.setResidueNumber(resNum)
        a.setChainID(chain)
        // make a copy of the atom to not re-parent the atom in the original assembly
        Atom aCopy = new Atom(a.xyzIndex, a, a.getXYZ(new double[3]), resNum, chain, Character.toString(chain))
        newAssembly.addMSNode(aCopy)
      }
      resNum++
    }

    // add the ions
    for (Atom atom : ions) {
      atom.setHetero(true)
      if (atom.getType() == naType) { // todo - could make list of all ion types (see above)
        atom.setResName("NA")
      } else if (atom.getType() == clType) {
        atom.setResName("CL")
      }
      atom.setResidueNumber(resNum)
      atom.setChainID(chain)
      // make a copy of the atom to not re-parent the atom in the original assembly
      Atom aCopy = new Atom(atom.xyzIndex, atom, atom.getXYZ(new double[3]), resNum, chain, Character.toString(chain))
      newAssembly.addMSNode(aCopy)
      resNum++
    }

    return newAssembly
  }

  // Creates a Polymer object for a DNA chain given an H5' to grow from
  private Polymer addDNAChain(Atom startAtom, Character chain) {
    // Set up new Polymer (chain) for a DNA strand
    Polymer polymer = new Polymer(chain, chain.toString())

    // atom to start building residue from
    Atom atom = startAtom

    int resNum = 1
    while (atom != null) {
      // atom list for current residue
      List<Atom> currResidue = new ArrayList<>()
      // add starting atom current residue atom list
      currResidue.add(atom)
      // atom list for next residue - will contain the next residue's starting atom (phosphorus)
      List<Atom> nextResidue = new ArrayList<>()
      // recursively find all atoms in the residue
      findDNAResidueAtoms(atom.index, atom, currResidue, nextResidue)
      // create a new residue
      Residue residue = new Residue(resNum, Residue.ResidueType.NA)
      residue.setChainID(chain)

      // add copied atoms to residue
      for (Atom a : currResidue) {
        int i = a.xyzIndex
        Atom aCopy = copiedAtomMap.get(i as Integer)
        aCopy.setResidueNumber(resNum)
        aCopy.setChainID(chain)
        aCopy.setSegID(chain as String)
        residue.addMSNode(aCopy)
      }

      String code = getDNAResidueCode(residue)
      residue.setName(code)

      // add residue to polymer
      polymer.addMSNode(residue)

      resNum++

      if (nextResidue.size() == 1) {
        atom = nextResidue.get(0) // start next residue at phosphorus
      } else {
        atom = null // end loop
      }
    }

    return polymer
  }

  // Creates a Polymer object for an amino acid chain given an NH3+ to grow from
  private Polymer addProteinChain(Atom startAtom, Character chain) {
    // Set up new Polymer (chain) for a protein strand
    Polymer polymer = new Polymer(chain, chain.toString())

    // atom to start building residue from
    Atom atom = startAtom

    int resNum = 1
    while (atom != null) {
      // atom list for current residue
      List<Atom> currResidue = new ArrayList<>()
      // add starting atom current residue atom list
      currResidue.add(atom)
      // atom list for next residue - will contain the next residue's starting atom (nitrogen)
      List<Atom> nextResidue = new ArrayList<>()
      // recursively find all atoms in the residue
      findAAResidueAtoms(atom.index, atom, currResidue, nextResidue)
      logger.info("NEXTRESSIZE: " + nextResidue.size() + " ; atoms: " + nextResidue)
      // create a new residue
      Residue residue = new Residue(resNum, Residue.ResidueType.AA)
      residue.setChainID(chain)

      // add copied atoms to residue
      for (Atom a : currResidue) {
        int i = a.xyzIndex
        Atom aCopy = copiedAtomMap.get(i as Integer)
        aCopy.setResidueNumber(resNum)
        aCopy.setChainID(chain)
        aCopy.setSegID(chain as String)
        residue.addMSNode(aCopy)
      }

      String code = getAAResidueCode(residue) // todo
      residue.setName(code)

      // add residue to polymer
      polymer.addMSNode(residue)

      resNum++

      if (nextResidue.size() == 1) {
        atom = nextResidue.get(0) // start next residue at peptide bond nitrogen
      } else {
        atom = null // end loop
      }
    }

    return polymer
  }

  // Recursive atom finding for an amino acid residue
  private void findAAResidueAtoms(int startingIndex, Atom atom, List<Atom> currResidue, List<Atom> nextResidue) {
    List<Atom> bondedAtoms = new ArrayList<>(atom.get12List())
    // If start atom (nitrogen that is part of the peptide bond
    if (atom.atomicNumber == 7 && isPeptideBond(atom)) {
      Atom cpep = null
      for (Atom a : bondedAtoms) {
        // look for the C atom in the peptide bond to remove it (part of the previous residue)
        if (a.atomicNumber == 6 && isPeptideBond(a)) {
          cpep = a
        }
      }
      if (cpep != null) {
        bondedAtoms.remove(cpep)
      } else {
        logger.severe("Missing C from the peptide bond.")
      }
    }
    for (Atom ba : bondedAtoms) {
      // Do not add the next amino acid's N of peptide bond to this residue. Add it to the next and stop recursion.
      if (!currResidue.contains(ba) && !(ba.atomicNumber == 7 && isPeptideBond(ba))) {
        currResidue.add(ba)
        findAAResidueAtoms(startingIndex, ba, currResidue, nextResidue)
      } else if ((ba.atomicNumber == 7 && isPeptideBond(ba))  && !nextResidue.contains(ba) && ba.index != startingIndex) {
        nextResidue.add(ba)
      }
    }
  }

  // Check to see if the atom is part of a peptide bond
  private static boolean isPeptideBond(Atom atom) {
    List<Atom> bondedAtoms = new ArrayList<>(atom.get12List())
    if (atom.atomicNumber == 7) {
      // The N atom of the bond should be bound to c-alpha (CA), hydrogen (HN), and peptide carbon (C)
      for (Atom a : bondedAtoms) {
        String atomName = a.atomType.name
//        String atomName = a.name
        if (atomName != "CA" && atomName != "HN" && atomName != "C") {
          return false
        }
      }
    } else if (atom.atomicNumber == 6) {
      // The C atom of the bond should be bound to c-alpha (CA), oxygen (O), and peptide nitrogen (N)
      for (Atom a : bondedAtoms) {
        String atomName = a.atomType.name
//        String atomName = a.name
        if (atomName != "CA" && atomName != "O" && atomName != "N") {
          return false
        }
      }
    } else {
      // A peptide bond only consists of nitrogen or carbon
      return false
    }

    // make sure there is only 3 bonds
    if (bondedAtoms.size() != 3) {
      return false
    }

    return true
  }

  // Recursive atom finding for a DNA residue
  private void findDNAResidueAtoms(int startingIndex, Atom atom, List<Atom> currResidue, List<Atom> nextResidue) {
    List<Atom> bondedAtoms = new ArrayList<>(atom.get12List())
    if (atom.atomicNumber == 15) { // 15 is phosphorus
      Atom o3 = null
      for (Atom a : bondedAtoms) {
        if (a.getType() == 338) { // 338 is Deoxyribose O3 - was part of previous residue
          o3 = a
        }
      }
      if (o3 != null) {
        bondedAtoms.remove(o3)
      } else {
        logger.severe("Missing O3 from Phosphorus bond list.")
      }
    }
    for (Atom ba : bondedAtoms) {
      if (!currResidue.contains(ba) && ba.atomicNumber != 15) { // 15 is phosphorus
        currResidue.add(ba)
        findDNAResidueAtoms(startingIndex, ba, currResidue, nextResidue)
      } else if (ba.atomicNumber == 15 && !nextResidue.contains(ba) && ba.index != startingIndex) {
        nextResidue.add(ba)
      }
    }
  }

  // Get amino acid residue name
  private String getAAResidueCode(Residue residue) {
    // Count the number of each element in the side chain and put into string (e.g., MET is "S1C3")
    List<Atom> atoms = new ArrayList<>(residue.getAtomList(true))
    logger.info("BEFORESIZE: " + atoms.size())
    logger.info("INITATOMS: " + atoms)
    List<String> backboneNames = ["N", "CA", "C", "O", "HN", "HA"] // backbone atom names, todo - ensure correct for proline
//    for (Atom atom : atoms) {
//      String atomName = atom.atomType.name
//      // Remove the atom from atom list if it is a backbone atom, only want sidechain
//      if (backboneNames.contains(atomName)) {
//        atoms.remove(atom)
//      }
//    }
//    atoms.removeIf {a -> backboneNames.contains(a.atomType.name)}
    atoms.removeIf {a -> backboneNames.contains(a.name)}
    logger.info("AFTERSIZE: " + atoms.size())

    // If removing the backbone atoms removes all atoms, the amino acid was glycine
    if (atoms.size() == 0) {
      return "GLY"
    }

    // Count S, O, N, and C atoms - excluding H from side chain uniqueness
    int nS = 0
    int nO = 0
    int nN = 0
    int nC = 0
    for (Atom atom : atoms) {
      if (atom.atomicNumber == 16) {
        nS++
      } else if (atom.atomicNumber == 8) {
        nO++
      } else if (atom.atomicNumber == 7) {
        nN++
      } else if (atom.atomicNumber == 6) {
        nC++
      }
    }

    StringBuilder sb = new StringBuilder()
    if (nS > 0) {
      sb.append("S" + nS)
    }
    if (nO > 0) {
      sb.append("O" + nO)
    }
    if (nN > 0) {
      sb.append("N" + nN)
    }
    if (nC > 0) {
      sb.append("C" + nC)
    }

    String code = AminoAcidUtils.sidechainStoichiometry.get(sb.toString())

    logger.info("ATOMS: " + atoms + " ; SS: " + sb + " ; CODE: " + code)

    if (code != null && code != "") {
      return code
    }

    // TODO:
    // Can't do Proline/Valine, and Isoleucine/leucine
    //    sidechainStoichiometry.put("C3", "1"); // Proline / Valine
    //    sidechainStoichiometry.put("C4", "2"); // (Iso)Leucine
    //    sidechainStoichiometry.put("O2N5C7", "3"); // DNA Gaunine / RNA Adenine

    logger.warning("UNMATCHED AMINO ACID: " + residue.toFormattedString(true, true))
    return ""
  }

  // Get residue name
  private String getDNAResidueCode(Residue residue) {
    List<Atom> atoms = residue.getAtomList(true)
    List<Integer> atomTypes = new ArrayList<>()
    for (Atom atom : atoms) {
      atomTypes.add(atom.getType())
    }
    atomTypes.sort()

    String[] dnaNames = ["Deoxythymidine", "Deoxyadenosine", "Deoxyguanosine", "Deoxycytidine"]
    for (String name : dnaNames) {
      // Work on a copy so the original dictionary entry is not mutated,
      // and remove by value (not index) explicitly.
      List<Integer> bioTypes = new ArrayList<>(moleculeAtomTypeDict.get(name))
      bioTypes.sort()

      // check if atom types match biotypes off the bat
      if (atomTypes == bioTypes) {
        return nameToCode(name)
      } else if (atomTypes.contains((Object) 342 as Integer)) {
        // if residue is a 5' end.. remove regular O5' and need to add H5T and O5' (terminal)
        bioTypes.remove((Object) 341) // remove regular O5'
        bioTypes.add((Object) 342 as Integer) // H5T - 5-Hydroxyl DNA
        bioTypes.add((Object) 348 as Integer) // O5' - 5-Hydroxyl DNA
        if (atomTypes == bioTypes) {
          return nameToCode(name)
        }
      } else if (atomTypes.contains((Object) 339 as Integer)) {
        // if residue is a 3' end.. remove regular O3' and need to add H3T and O3' (terminal) along with phosphate
        bioTypes.remove((Object) 338) // remove regular O3'
        bioTypes.add((Object) 343 as Integer) // P - Phosphodiester DNA
        bioTypes.add((Object) 344 as Integer) // OP - Phosphodiester DNA
        bioTypes.add((Object) 344 as Integer) // OP - Phosphodiester DNA
        bioTypes.add((Object) 339 as Integer) // H3T - 3-Hydroxyl DNA
        bioTypes.add((Object) 347 as Integer) // O3' - 3-Hydroxyl DNA
        bioTypes.sort()
        if (atomTypes == bioTypes) {
          return nameToCode(name)
        }
      } else {
        // not 5' or 3' terminal residue, so just add phosphate group
        bioTypes.add((Object) 343 as Integer) // P - Phosphodiester DNA
        bioTypes.add((Object) 344 as Integer) // OP - Phosphodiester DNA
        bioTypes.add((Object) 344 as Integer) // OP - Phosphodiester DNA
        if (atomTypes == bioTypes) {
          return nameToCode(name)
        }
      }

      // final check - need to add the other 2 H7s if the residue is deoxythymidine (bioTypes contains 1 H7, need 3)
      if (name == "Deoxythymidine" && atomTypes.contains((Object) 295 as Integer)) {
        bioTypes.add((Object) 295 as Integer) // H7 - Deoxythymidine
        bioTypes.add((Object) 295 as Integer) // H7 - Deoxythymidine
        bioTypes.sort()
        if (atomTypes == bioTypes) {
          return nameToCode(name)
        }
      }
    }

    // todo - no support for a single nucleotide.. won't handle both 5'-OH and  3'-OH
    // todo - could use stoichiometry like amino acids - not for now because it works

    logger.warning("UNMATCHED RESIDUE: " + residue.toFormattedString(true, true))
    return ""
  }

  private static String nameToCode(String name) {
    String code
    switch (name) {
      case "Deoxythymidine" -> {
        code = "DTY"
      }
      case "Deoxyadenosine" -> {
        code = "DAD"
      }
      case "Deoxyguanosine"-> {
        code = "DGU"
      }
      case "Deoxycytidine" -> {
        code = "DCY"
      }
      default -> {
        code = "Unknown"
      }
    }
    return code
  }

}
