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
package ffx.potential;

import ffx.crystal.Crystal;
import ffx.crystal.SpaceGroup;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly.FractionalMode;
import ffx.potential.bonded.Atom;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.toDegrees;

/**
 * This class computes the energy and Cartesian coordinate gradient, plus finite difference
 * derivatives of lattice parameters.
 *
 * @author Jooyeon Park
 * @since 1.0
 */
public class XtalEnergy implements Potential {

  /**
   * MolecularAssembly to compute the crystal energy for.
   */
  private final MolecularAssembly molecularAssembly;
  /**
   * The ForceFieldEnergy to use.
   */
  private final ForceFieldEnergy forceFieldEnergy;
  /**
   * Array of active atoms.
   */
  private final Atom[] activeAtoms;
  /**
   * Number of active atoms.
   */
  private final int nActive;
  /**
   * Atomic coordinates
   */
  private final double[] xyz;
  private final double[] gr;
  private final int nParams;

  private final Crystal crystal;
  private final boolean latticeOnly;
  private final VARIABLE_TYPE[] type;
  private final double[] mass;
  private final Crystal unitCell;
  private double[] scaling;
  private double totalEnergy;
  private FractionalMode fractionalMode = FractionalMode.OFF;

  /**
   * Constructor for XtalEnergy.
   *
   * @param forceFieldEnergy  a {@link ffx.potential.ForceFieldEnergy} object.
   * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
   */
  public XtalEnergy(ForceFieldEnergy forceFieldEnergy, MolecularAssembly molecularAssembly) {
    this(forceFieldEnergy, molecularAssembly, false);
  }

  /**
   * Constructor for XtalEnergy.
   *
   * @param forceFieldEnergy  a {@link ffx.potential.ForceFieldEnergy} object.
   * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
   * @param latticeOnly       if true, only lattice parameters are optimized.
   */
  public XtalEnergy(ForceFieldEnergy forceFieldEnergy, MolecularAssembly molecularAssembly, boolean latticeOnly) {
    this.forceFieldEnergy = forceFieldEnergy;
    this.molecularAssembly = molecularAssembly;
    this.latticeOnly = latticeOnly;

    if (!latticeOnly) {
      Atom[] atoms = molecularAssembly.getAtomArray();
      List<Atom> active = new ArrayList<>();
      for (Atom a : atoms) {
        if (a.isActive()) {
          active.add(a);
        }
      }
      nActive = active.size();
      activeAtoms = active.toArray(new Atom[0]);
      xyz = new double[3 * nActive];
      gr = new double[3 * nActive];
    } else {
      nActive = 0;
      activeAtoms = null;
      xyz = null;
      gr = null;
    }

    nParams = 3 * nActive + 6;
    crystal = forceFieldEnergy.getCrystal();
    unitCell = crystal.getUnitCell();
    type = new VARIABLE_TYPE[nParams];
    mass = new double[nParams];

    // Load atomic masses and variable types for atomic coordinates.
    int index = 0;
    for (int i = 0; i < nActive; i++) {
      double m = activeAtoms[i].getMass();
      mass[index] = m;
      mass[index + 1] = m;
      mass[index + 2] = m;
      type[index] = VARIABLE_TYPE.X;
      type[index + 1] = VARIABLE_TYPE.Y;
      type[index + 2] = VARIABLE_TYPE.Z;
      index += 3;
    }
    // Load variable types for lattice parameters.
    for (int i = nActive * 3; i < nActive * 3 + 6; i++) {
      mass[i] = 1.0;
      type[i] = VARIABLE_TYPE.OTHER;
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean destroy() {
    return forceFieldEnergy.destroy();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energy(double[] x) {

    // Un-scale coordinates if applicable.
    unscaleCoordinates(x);

    // Set atomic coordinates & lattice parameters.
    setCoordinates(x);

    // Compute energy.
    totalEnergy = forceFieldEnergy.energy(false, false);

    // Scale coordinates if applicable.
    scaleCoordinates(x);

    return totalEnergy;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energyAndGradient(double[] x, double[] g) {
    // Un-scale coordinates if applicable.
    unscaleCoordinates(x);

    // Set atomic coordinates & lattice parameters.
    setCoordinates(x);

    if (latticeOnly) {
      // Only lattice parameters are optimized.
      totalEnergy = forceFieldEnergy.energy(false, false);
    } else {
      // Compute energy and gradient.
      totalEnergy = forceFieldEnergy.energyAndGradient(xyz, gr);

    }

    // Both coordinates and gradient are scaled if applicable.
    packGradient(x, g);

    // Calculate finite-difference partial derivatives of lattice parameters.
    unitCellParameterDerivatives(x, g);

    return totalEnergy;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getAcceleration(double[] acceleration) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getCoordinates(double[] x) {
    int n = getNumberOfVariables();
    if (x == null || x.length < n) {
      x = new double[n];
    }
    int index = 0;
    for (int i = 0; i < nActive; i++) {
      Atom a = activeAtoms[i];
      x[index] = a.getX();
      index++;
      x[index] = a.getY();
      index++;
      x[index] = a.getZ();
      index++;
    }
    x[index] = unitCell.a;
    index++;
    x[index] = unitCell.b;
    index++;
    x[index] = unitCell.c;
    index++;
    x[index] = unitCell.alpha;
    index++;
    x[index] = unitCell.beta;
    index++;
    x[index] = unitCell.gamma;
    return x;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public STATE getEnergyTermState() {
    return forceFieldEnergy.getEnergyTermState();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setEnergyTermState(STATE state) {
    forceFieldEnergy.setEnergyTermState(state);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getMass() {
    return mass;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfVariables() {
    return nParams;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getPreviousAcceleration(double[] previousAcceleration) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getScaling() {
    return scaling;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setScaling(double[] scaling) {
    this.scaling = scaling;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getTotalEnergy() {
    return totalEnergy;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public VARIABLE_TYPE[] getVariableTypes() {
    return type;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getVelocity(double[] velocity) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setAcceleration(double[] acceleration) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /**
   * setFractionalCoordinateMode.
   *
   * @param fractionalMode a {@link ffx.potential.MolecularAssembly.FractionalMode} object.
   */
  public void setFractionalCoordinateMode(FractionalMode fractionalMode) {
    this.fractionalMode = fractionalMode;
    molecularAssembly.setFractionalMode(fractionalMode);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setPreviousAcceleration(double[] previousAcceleration) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setVelocity(double[] velocity) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  /**
   * Use finite-differences to compute unit cell derivatives.
   *
   * @param x Coordinates and unit cell parameters.
   * @param g gradient.
   */
  private void unitCellParameterDerivatives(double[] x, double[] g) {
    int index = 3 * nActive;
      double eps = 1.0e-5;
      double deps = toDegrees(eps);
      SpaceGroup spaceGroup = crystal.spaceGroup;
      switch (spaceGroup.latticeSystem) {
        case TRICLINIC_LATTICE -> {
          g[index] = finiteDifference(x, index, eps);
          index++;
          g[index] = finiteDifference(x, index, eps);
          index++;
          g[index] = finiteDifference(x, index, eps);
          index++;
          g[index] = finiteDifference(x, index, deps);
          index++;
          g[index] = finiteDifference(x, index, deps);
          index++;
          g[index] = finiteDifference(x, index, deps);
        }
        case MONOCLINIC_LATTICE -> {
          // alpha = gamma = 90
          g[index] = finiteDifference(x, index, eps);
          index++;
          g[index] = finiteDifference(x, index, eps);
          index++;
          g[index] = finiteDifference(x, index, eps);
          index++;
          g[index] = 0.0;
          index++;
          g[index] = finiteDifference(x, index, deps);
          index++;
          g[index] = 0.0;
        }
        case ORTHORHOMBIC_LATTICE -> {
          // alpha = beta = gamma = 90
          g[index] = finiteDifference(x, index, eps);
          index++;
          g[index] = finiteDifference(x, index, eps);
          index++;
          g[index] = finiteDifference(x, index, eps);
          index++;
          g[index] = 0.0;
          index++;
          g[index] = 0.0;
          index++;
          g[index] = 0.0;
        }
        case TETRAGONAL_LATTICE, // a = b, alpha = beta = gamma = 90
             HEXAGONAL_LATTICE // a = b, alpha = beta = 90, gamma = 120
            -> {
          g[index] = finiteDifference2(x, index, index + 1, eps);
          index++;
          g[index] = g[index - 1];
          index++;
          g[index] = finiteDifference(x, index, eps);
          index++;
          g[index] = 0.0;
          index++;
          g[index] = 0.0;
          index++;
          g[index] = 0.0;
        }
        case RHOMBOHEDRAL_LATTICE -> {
          // a = b = c, alpha = beta = gamma
          g[index] = finiteDifference3(x, index, index + 1, index + 2, eps);
          index++;
          g[index] = g[index - 1];
          index++;
          g[index] = g[index - 2];
          index++;
          g[index] = finiteDifference3(x, index, index + 1, index + 2, deps);
          index++;
          g[index] = g[index - 1];
          index++;
          g[index] = g[index - 2];
        }
        // a = b, alpha = beta = 90, gamma = 120
        case CUBIC_LATTICE -> {
          // a = b = c, alpha = beta = gamma = 90
          g[index] = finiteDifference3(x, index, index + 1, index + 2, eps);
          index++;
          g[index] = g[index - 1];
          index++;
          g[index] = g[index - 2];
          index++;
          g[index] = 0.0;
          index++;
          g[index] = 0.0;
          index++;
          g[index] = 0.0;
        }
      }

    // Scale finite-difference partial derivatives of lattice parameters.
    // Coordinate scaling has already been applied.
    if (scaling != null) {
      index = 3 * nActive;
      for (int i = 0; i < 6; i++) {
        g[index] /= scaling[index];
        index++;
      }
    }
  }

  /**
   * Calculate finite-difference derivative for any parameter.
   *
   * @param x     Coordinates and unit cell parameters.
   * @param index Parameter index.
   * @param eps   Step size.
   * @return Finite-difference derivative.
   */
  private double finiteDifference(double[] x, int index, double eps) {
    double scale = 1.0;
    if (scaling != null) {
      scale = scaling[index];
    }
    double x1 = x[index];
    final double param = x1 / scale;

    x[index] = (param + eps / 2.0) * scale;
    double ePlus = energy(x);
    x[index] = (param - eps / 2.0) * scale;
    double eMinus = energy(x);

    x[index] = x1;

    return (ePlus - eMinus) / eps;
  }

  /**
   * @param x      Coordinates and unit cell parameters.
   * @param index1 Parameter index 1.
   * @param index2 Parameter index 2.
   * @param eps    Step size.
   * @return Finite-difference derivative.
   */
  private double finiteDifference2(double[] x, int index1, int index2, double eps) {
    double scale1 = 1.0;
    double scale2 = 1.0;

    if (scaling != null) {
      scale1 = scaling[index1];
      scale2 = scaling[index2];
    }

    final double x1 = x[index1];
    final double x2 = x[index2];
    final double param1 = x1 / scale1;
    final double param2 = x2 / scale2;
    x[index1] = (param1 + eps / 2.0) * scale1;
    x[index2] = (param2 + eps / 2.0) * scale2;
    double ePlus = energy(x);
    x[index1] = (param1 - eps / 2.0) * scale1;
    x[index2] = (param2 - eps / 2.0) * scale2;
    double eMinus = energy(x);

    x[index1] = x1;
    x[index2] = x2;

    return (ePlus - eMinus) / eps;
  }

  /**
   * @param x      Coordinates and unit cell parameters.
   * @param index1 Parameter index 1.
   * @param index2 Parameter index 2.
   * @param index3 Parameter index 3.
   * @param eps    Step size.
   * @return finite-difference derivative.
   */
  private double finiteDifference3(double[] x, int index1, int index2, int index3, double eps) {
    double scale1 = 1.0;
    double scale2 = 1.0;
    double scale3 = 1.0;
    if (scaling != null) {
      scale1 = scaling[index1];
      scale2 = scaling[index2];
      scale3 = scaling[index3];
    }

    final double x1 = x[index1];
    final double x2 = x[index2];
    final double x3 = x[index3];

    final double param1 = x1 / scale1;
    final double param2 = x2 / scale2;
    final double param3 = x3 / scale3;
    x[index1] = (param1 + eps / 2.0) * scale1;
    x[index2] = (param2 + eps / 2.0) * scale2;
    x[index3] = (param3 + eps / 2.0) * scale3;
    final double ePlus = energy(x);
    x[index1] = (param1 - eps / 2.0) * scale1;
    x[index2] = (param2 - eps / 2.0) * scale2;
    x[index3] = (param3 - eps / 2.0) * scale3;
    final double eMinus = energy(x);

    x[index1] = x1;
    x[index2] = x2;
    x[index3] = x3;

    return (ePlus - eMinus) / eps;
  }

  /**
   * Apply scaling for the optimizer if applicable.
   *
   * @param x Coordinates and unit cell parameters.
   * @param g Atomic coordinate gradient.
   */
  private void packGradient(double[] x, double[] g) {
    // Scale fractional coordinates and gradient.
    if (scaling != null) {
      int len = x.length;
      for (int i = 0; i < len; i++) {
        if (fractionalMode == FractionalMode.OFF) {
          g[i] /= scaling[i];
        } else {
          // If we're maintaining fractional coordinates,
          // the atomic coordinate derivatives can be zero.
          g[i] = 0.0;
        }
        x[i] *= scaling[i];
      }
    }
  }

  /**
   * Sets atomic coordinates and lattice parameters.
   *
   * @param x First 3*nActive parameters are coordinates, next 6 are x parameters.
   */
  private void setCoordinates(double[] x) {
    assert (x != null);

    // Before applying new lattice parameters, store factional coordinates.
    if (fractionalMode != FractionalMode.OFF) {
      molecularAssembly.computeFractionalCoordinates();
    }

    int index = nActive * 3;
    double a = x[index];
    double b = x[index + 1];
    double c = x[index + 2];
    double alpha = x[index + 3];
    double beta = x[index + 4];
    double gamma = x[index + 5];

    SpaceGroup spaceGroup = crystal.spaceGroup;

    // Enforce the lattice system.
    switch (spaceGroup.latticeSystem) {
      case TRICLINIC_LATTICE -> {
      }
      case MONOCLINIC_LATTICE -> {
        // alpha = gamma = 90
        alpha = 90.0;
        gamma = 90.0;
      }
      case ORTHORHOMBIC_LATTICE -> {
        // alpha = beta = gamma = 90
        alpha = 90.0;
        beta = 90.0;
        gamma = 90.0;
      }
      case TETRAGONAL_LATTICE -> {
        // a = b, alpha = beta = gamma = 90
        double ab = 0.5 * (a + b);
        a = ab;
        b = ab;
        alpha = 90.0;
        beta = 90.0;
        gamma = 90.0;
      }
      case RHOMBOHEDRAL_LATTICE -> {
        // a = b = c, alpha = beta = gamma.
        double abc = (a + b + c) / 3.0;
        a = abc;
        b = abc;
        c = abc;
        double angles = (alpha + beta + gamma) / 3.0;
        alpha = angles;
        beta = angles;
        gamma = angles;
      }
      case HEXAGONAL_LATTICE -> {
        // a = b, alpha = beta = 90 && gamma = 120
        double ab = 0.5 * (a + b);
        a = ab;
        b = ab;
        alpha = 90.0;
        beta = 90.0;
        gamma = 120.0;
      }
      case CUBIC_LATTICE -> {
        // a = b = c, alpha = beta = gamma = 90
        double abc = (a + b + c) / 3.0;
        a = abc;
        b = abc;
        c = abc;
        alpha = 90.0;
        beta = 90.0;
        gamma = 90.0;
      }
    }

    crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma);
    forceFieldEnergy.setCrystal(crystal);

    // Use the atomic coordinates from the optimizer.
    if (fractionalMode == FractionalMode.OFF) {
      index = 0;
      for (int i = 0; i < nActive; i++) {
        Atom atom = activeAtoms[i];
        double xx = x[index];
        double yy = x[index + 1];
        double zz = x[index + 2];
        xyz[index] = xx;
        xyz[index + 1] = yy;
        xyz[index + 2] = zz;
        index += 3;
        atom.moveTo(xx, yy, zz);
      }
    } else {
      // Maintain fractional coordinates.
      molecularAssembly.moveToFractionalCoordinates();
      index = 0;
      for (int i = 0; i < nActive; i++) {
        Atom atom = activeAtoms[i];
        double xx = atom.getX();
        double yy = atom.getY();
        double zz = atom.getZ();
        xyz[index] = xx;
        xyz[index + 1] = yy;
        xyz[index + 2] = zz;
        index += 3;
      }
    }
  }
}
