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
package ffx.ui.behaviors;

import java.awt.event.MouseEvent;
import org.jogamp.java3d.Bounds;
import org.jogamp.java3d.BranchGroup;
import org.jogamp.java3d.Canvas3D;
import org.jogamp.java3d.Transform3D;
import org.jogamp.java3d.TransformGroup;
import org.jogamp.java3d.utils.picking.PickResult;

/**
 * The PickTranslateBehavior class implements a translation behavior on a picked scenegraph object.
 *
 * @author Michael J. Schnieders
 */
public class PickTranslateBehavior extends PickMouseBehavior implements MouseBehaviorCallback {

  public MouseTranslate translate;
  private PickingCallback callback = null;
  private TransformGroup currentTG;

  /**
   * Constructor for PickTranslateBehavior.
   *
   * @param root a {@link org.jogamp.java3d.BranchGroup} object.
   * @param canvas a {@link org.jogamp.java3d.Canvas3D} object.
   * @param bounds a {@link org.jogamp.java3d.Bounds} object.
   * @param VPTG a {@link org.jogamp.java3d.TransformGroup} object.
   * @param pickMode a int.
   */
  public PickTranslateBehavior(
      BranchGroup root, Canvas3D canvas, Bounds bounds, TransformGroup VPTG, int pickMode) {
    super(canvas, root, bounds);
    translate = new MouseTranslate(MouseBehavior.MANUAL_WAKEUP, VPTG);
    translate.setTransformGroup(currGrp);
    translate.setFactor(0.1);
    currGrp.addChild(translate);
    translate.setSchedulingBounds(bounds);
    setSchedulingBounds(bounds);
    pickCanvas.setMode(pickMode);
  }

  /**
   * Return the pickMode component of this PickTranslateBehavior.
   *
   * @return a int.
   */
  public int getPickMode() {
    return pickCanvas.getMode();
  }

  /**
   * Sets the pickMode component of this PickTranslateBehavior to the value of the passed
   * pickMode. @param pickMode the pickMode to be copied.
   *
   * @param pickMode a int.
   */
  public void setPickMode(int pickMode) {
    pickCanvas.setMode(pickMode);
  }

  /**
   * Register the class @param callback to be called each time the picked object moves
   *
   * @param callback a {@link ffx.ui.behaviors.PickingCallback} object.
   */
  public void setupCallback(PickingCallback callback) {
    this.callback = callback;
    if (callback == null) {
      translate.setupCallback(null);
    } else {
      translate.setupCallback(this);
    }
  }

  /** {@inheritDoc} */
  public void transformChanged(int type, Transform3D transform) {
    callback.transformChanged(PickingCallback.TRANSLATE, currentTG);
  }

  /** {@inheritDoc} */
  public void transformClicked(int type, Transform3D transform) {
    callback.transformClicked(PickingCallback.TRANSLATE, currentTG);
  }

  /** {@inheritDoc} */
  public void transformDoubleClicked(int type, Transform3D transform) {
    callback.transformDoubleClicked(PickingCallback.TRANSLATE, currentTG);
  }

  /**
   * Update the scene to manipulate any nodes. This is not meant to be called by users. Behavior
   * automatically calls this. You can call this only if you know what you are doing.
   *
   * @param xpos Current mouse X pos.
   * @param ypos Current mouse Y pos.
   */
  public void updateScene(int xpos, int ypos) {
    if ((mevent.getModifiersEx() & MouseEvent.BUTTON3_DOWN_MASK) == MouseEvent.BUTTON3_DOWN_MASK) {
      pickCanvas.setShapeLocation(xpos, ypos);
      PickResult r = pickCanvas.pickClosest();
      if (r != null) {
        if (callback != null) {
          callback.transformChanged(PickingCallback.NO_PICK, null);
        }
      }
    }
  }
}
