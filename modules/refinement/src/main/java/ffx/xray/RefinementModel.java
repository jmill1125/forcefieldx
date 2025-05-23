// ******************************************************************************
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
// ******************************************************************************
package ffx.xray;

import static java.lang.String.format;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * RefinementModel class.
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public class RefinementModel {

  private static final Logger logger = Logger.getLogger(RefinementModel.class.getName());
  /** An atom list. */
  private final List<Atom> totalAtomList;
  /** An atom array. */
  private final Atom[] totalAtomArray;
  /** An atom list. */
  private final List<Atom> activeAtomList;
  /** An atom array. */
  private final Atom[] activeAtomArray;
  /** Map between atom in different molecular assemblies. */
  private final List<List<Integer>> xIndex;

  private final List<List<Residue>> altResidues;
  private final List<List<Molecule>> altMolecules;

  /**
   * Constructor for RefinementModel.
   *
   * @param molecularAssemblies an array of {@link ffx.potential.MolecularAssembly} objects.
   */
  public RefinementModel(MolecularAssembly[] molecularAssemblies) {
    this(molecularAssemblies, false);
  }

  /**
   * Constructor for RefinementModel.
   *
   * @param molecularAssemblies an array of {@link ffx.potential.MolecularAssembly} objects.
   * @param refineMolOcc a boolean.
   */
  @SuppressWarnings("unchecked")
  public RefinementModel(MolecularAssembly[] molecularAssemblies, boolean refineMolOcc) {
    List<Atom> atomList;

    // Build alternate conformer list for occupancy refinement (if necessary).
    altResidues = new ArrayList<>();
    altMolecules = new ArrayList<>();
    List<MSNode> nodeList0 = molecularAssemblies[0].getNodeList();
    List<Residue> tempResidues = null;
    List<Molecule> tempMolecules = null;
    boolean alternateConformer;

    // By residue/molecule.
    for (int i = 0; i < nodeList0.size(); i++) {
      alternateConformer = false;
      MSNode iNode = nodeList0.get(i);

      // First set up alternate residue restraint list.
      for (Atom a : iNode.getAtomList()) {
        if (!a.getUse()) {
          continue;
        }
        if (a.getAltLoc() == null) {
          logger.severe(
              format(
                  " Atom %s has a null alternate location. Likely cause: attempting X-ray refinement from a .xyz file",
                  a));
        }
        if (!a.getAltLoc().equals(' ') || a.getOccupancy() < 1.0) {
          if (iNode instanceof Residue) {
            tempResidues = new ArrayList<>();
            tempResidues.add((Residue) iNode);
            altResidues.add(tempResidues);
            alternateConformer = true;
            break;
          } else if (iNode instanceof Molecule) {
            if (refineMolOcc) {
              tempMolecules = new ArrayList<>();
              tempMolecules.add((Molecule) iNode);
              altMolecules.add(tempMolecules);
            }
            alternateConformer = true;
            break;
          }
        }
      }
      if (alternateConformer) {
        for (int j = 1; j < molecularAssemblies.length; j++) {
          List<MSNode> nlist = molecularAssemblies[j].getNodeList();
          MSNode node = nlist.get(i);
          for (Atom a : node.getAtomList()) {
            if (!a.getUse()) {
              continue;
            }
            if (!a.getAltLoc().equals(' ') && !a.getAltLoc().equals('A')) {
              if (node instanceof Residue) {
                if (tempResidues != null) {
                  tempResidues.add((Residue) node);
                }
                break;
              } else if (node instanceof Molecule) {
                if (tempMolecules != null) {
                  tempMolecules.add((Molecule) node);
                }
                break;
              }
            }
          }
        }
      }
    }

    // For mapping between atoms between different molecular assemblies.
    // xIndex = new List[molecularAssemblies.length];
    xIndex = new ArrayList<>(molecularAssemblies.length);
    for (int i = 0; i < molecularAssemblies.length; i++) {
      // xIndex[i] = new ArrayList<>();
      xIndex.add(new ArrayList<>());
    }

    // Create the atomList to be used for SF calculations.
    totalAtomList = new ArrayList<>();

    // Create the activeAtomList (i.e. atoms that can move).
    activeAtomList = new ArrayList<>();

    // Root list.
    atomList = molecularAssemblies[0].getAtomList();
    int index = 0;
    for (Atom a : atomList) {
      if (!a.getUse()) {
        continue;
      }
      totalAtomList.add(a);
//      a.setFormFactorIndex(index);
//      xIndex[0].add(index);
//      index++;
      if (a.isActive()) {
        activeAtomList.add(a);
        a.setFormFactorIndex(index);
        // xIndex[0].add(index);
        xIndex.get(0).add(index);
        index++;
      }
    }

    // Now add cross-references to root and any alternate atoms not in root
    for (int i = 1; i < molecularAssemblies.length; i++) {
      atomList = molecularAssemblies[i].getAtomList();
      for (Atom a : atomList) {
        if (!a.getUse()) {
          continue;
        }
        Atom root = molecularAssemblies[0].findAtom(a);
        if (root != null && root.getAltLoc().equals(a.getAltLoc())) {
//          xIndex[i].add(root.getFormFactorIndex());
//          a.setFormFactorIndex(root.getFormFactorIndex());
          if (a.isActive()) {
            // xIndex[i].add(root.getFormFactorIndex());
            xIndex.get(i).add(root.getFormFactorIndex());
            a.setFormFactorIndex(root.getFormFactorIndex());
          }
        } else {
//          xIndex[i].add(index);
//          index++;
          totalAtomList.add(a);
          if (a.isActive()) {
            activeAtomList.add(a);
            // xIndex[i].add(index);
            xIndex.get(i).add(index);
            index++;
          }
        }
      }
    }

    totalAtomArray = totalAtomList.toArray(new Atom[0]);
    activeAtomArray = activeAtomList.toArray(new Atom[0]);

    /*
     Make sure the occupancy values make sense, otherwise print warnings
     (since this could destabilize the refinement, should we error out?)
    */
    for (List<Residue> list : altResidues) {
      double tocc = 0.0;
      for (Residue r : list) {
        for (Atom a : r.getAtomList()) {
          if (a.getOccupancy() < 1.0 || a.getOccupancy() > 1.0) {
            tocc += a.getOccupancy();
            break;
          }
        }
      }
      if (list.size() == 1) {
        Residue r = list.get(0);
        logger.log(Level.INFO,
            "  Residue {0} is a single conformer with non-unity occupancy.\n  Occupancy will be refined independently.\n",
            r.getChainID() + " " + r);
      } else if (tocc < 1.0 || tocc > 1.0) {
        Residue r = list.get(0);
        logger.log(Level.INFO,
            "  Residue {0} occupancy does not sum to 1.0.\n  This should be fixed or checked due to possible instability in refinement.\n",
            r.getChainID() + " " + r);
      }
    }
  }

  /**
   * Getter for the field <code>activeAtomArray</code>.
   *
   * @return the activeAtomArray
   */
  public Atom[] getActiveAtomArray() {
    return activeAtomArray;
  }

  /**
   * Getter for the field <code>activeAtomList</code>.
   *
   * @return the activeAtomList
   */
  public List<Atom> getActiveAtomList() {
    return activeAtomList;
  }

  /**
   * Getter for the field <code>altMolecules</code>.
   *
   * @return the altMolecules
   */
  public List<List<Molecule>> getAltMolecules() {
    return altMolecules;
  }

  /**
   * Getter for the field <code>altResidues</code>.
   *
   * @return the altResidues
   */
  public List<List<Residue>> getAltResidues() {
    return altResidues;
  }

  /**
   * Getter for the field <code>totalAtomArray</code>.
   *
   * @return the totalAtomArray
   */
  public Atom[] getTotalAtomArray() {
    return totalAtomArray;
  }

  /**
   * Getter for the field <code>xIndex</code>.
   *
   * @return the xIndex
   */
  List<List<Integer>> getxIndex() {
    return xIndex;
  }

  /**
   * Getter for the field <code>totalAtomList</code>.
   *
   * @return the totalAtomList
   */
  List<Atom> getTotalAtomList() {
    return totalAtomList;
  }
}
