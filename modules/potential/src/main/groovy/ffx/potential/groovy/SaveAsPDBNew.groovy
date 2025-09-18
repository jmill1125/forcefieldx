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
import ffx.potential.bonded.Atom
import ffx.potential.bonded.MSNode
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.cli.PotentialScript
import ffx.potential.cli.SaveOptions
import ffx.potential.extended.ExtendedSystem
import ffx.potential.parameters.BioType
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XPHFilter
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
   * --esv Handle an extended system at the bottom of XYZ files using XPHFilter.
   */
  @Option(names = ['--esv'], paramLabel = "file", defaultValue = "",
      description = 'PDB file to build extended system from.')
  private String extended = ""

  /**
   * The final argument is an XYZ or ARC coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in XYZ or ARC format.')
  private String filename = null

  private Map<String, List<Integer>> moleculeAtomTypeDict = new HashMap<>()

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

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath()

    // CREATE BIOTYPE MOLECULE MAP - { (MOL1, AtomType[]), (MOL2, AtomType[]), ... }
    Map<String, BioType> bioTypes = activeAssembly.getForceField().getBioTypeMap()
    Map<String, ArrayList<BioType>> moleculeDict = new HashMap<>() // todo can delete but careful with logic
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
    int fivePrimeHType = moleculeAtomTypeDict.get("5-Hydroxyl DNA").get(1) // 0 = O5* ; 1 = H5T
    int waterType = moleculeAtomTypeDict.get("Water").get(0) // 0 = O ; 1 = H
    int naType = moleculeAtomTypeDict.get("Sodium Ion").get(0)
    int clType = moleculeAtomTypeDict.get("Chloride Ion").get(0)
    List<Atom> fivePrimeOs = new ArrayList<>()
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
        fivePrimeOs.add(atom)
      }
    }

    // Make all DNA chains
    // todo - assuming all DNA strands start with HO5' and end with HO3'
    int num = 0
    MSNode polymers = new MSNode()
    for (Atom atom : fivePrimeOs) {
      Character chain = 'X'
      if (num == 0) {
        chain = 'A'
      } else if (num == 1) {
        chain = 'B'
      } else {
        logger.severe("MORE THAN 2 CHAINS. CAN'T HANDLE YET.") // todo
      }
      Polymer polymer = addDNAChain(atom, chain) // todo - could also do for protein (and RNA)
      polymers.add(polymer)
      num++
    }

    // create new molecular assembly with polymers created above
    MolecularAssembly newAssembly = new MolecularAssembly("NewAssembly", polymers, activeAssembly.getForceField().getProperties())
    newAssembly.setForceField(activeAssembly.getForceField())

    // add the water molecules created above
    Character chain = 'C'
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
        newAssembly.addMSNode(a)
      }
      resNum++
    }

    // add the ions
    for (Atom atom : ions) {
      atom.setHetero(true)
      if (atom.getType() == naType) {
        atom.setResName("NA")
      } else if (atom.getType() == clType) {
        atom.setResName("CL")
      }
      atom.setResidueNumber(resNum)
      atom.setChainID(chain)
      newAssembly.addMSNode(atom)
      resNum++
    }

    String name = removeExtension(getName(filename)) + ".pdb"
    File saveFile = new File(name)

    saveFile = potentialFunctions.versionFile(saveFile)
    PDBFilter pdbFilter = new PDBFilter(saveFile, newAssembly, newAssembly.getForceField(), newAssembly.getProperties())

    saveFile.append(activeAssembly.getCrystal().toCRYST1()) // append header CRYST1 line
    pdbFilter.writeFile(saveFile, true)

    return this
  }

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

      // add atoms to residue
      for (Atom a : currResidue) {
        residue.addMSNode(a)
      }

      String code = getResidueCode(residue)
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

  // Attempt at Recursive atom finding
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

  // Get residue name
  private String getResidueCode(Residue residue) {
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

      if (atomTypes == bioTypes) {
        return nameToCode(name)
      } else if (atomTypes.contains((Object) 342 as Integer)) {
        bioTypes.remove((Object) 341) // remove regular O5' // todo - would it work for all four bases ??
        bioTypes.add((Object) 342 as Integer) // H5T - 5-Hydroxyl DNA
        bioTypes.add((Object) 348 as Integer) // O5' - 5-Hydroxyl DNA
        if (atomTypes == bioTypes) {
          return nameToCode(name)
        }
      } else if (atomTypes.contains((Object) 339 as Integer)) {
        // todo - no support for a single nucleotide.. won't handle both 5'-OH and  3'-OH

        bioTypes.remove((Object) 338) // remove regular O3' // todo - would it work for all four bases ??
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
        bioTypes.add((Object) 343 as Integer) // P - Phosphodiester DNA
        bioTypes.add((Object) 344 as Integer) // OP - Phosphodiester DNA
        bioTypes.add((Object) 344 as Integer) // OP - Phosphodiester DNA
        if (atomTypes == bioTypes) {
          return nameToCode(name)
        } else if (name == "Deoxythymidine" && atomTypes.contains((Object) 295 as Integer)) {
          bioTypes.add((Object) 295 as Integer) // H7 - Deoxythymidine
          bioTypes.add((Object) 295 as Integer) // H7 - Deoxythymidine
          bioTypes.sort()
          if (atomTypes == bioTypes) {
            return nameToCode(name)
          }
        }
      }
    }

    logger.info("UNMATCHED RESIDUE: " + residue.toFormattedString(true, true))
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
