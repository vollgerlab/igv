package org.igv.alignment.ma;

import org.igv.alignment.AlignmentCounts;
import org.igv.track.RenderContext;

import java.awt.*;


/**
 * Draws MA (molecular annotation) composition into the coverage bar at a single position, the
 * analog of {@link org.igv.alignment.mods.BaseModificationCoverageRenderer}.
 * <p>
 * Each type is drawn as its own bar rising from the bottom, in the fixed back to front order
 * {@link MaCounts#getTypesInDrawOrder()} defines, so FIRE ends up along the bottom.  Overlay rather
 * than stacking, because MA types are frequently nested rather than disjoint -- a FIRE is an MSP
 * that cleared an FDR threshold, so stacking them would report more than 100% of the reads at that
 * position.  Each bar's height is that type's own fraction of the reads covering the position,
 * which reads correctly for nested and disjoint types alike.
 */
public class MaCoverageRenderer {

    public static void drawAnnotations(RenderContext context,
                                       int pX,
                                       int pBottom,
                                       int dX,
                                       int barHeight,
                                       int pos,
                                       AlignmentCounts alignmentCounts) {

        MaCounts maCounts = alignmentCounts.getMaCounts();
        if (maCounts == null || maCounts.isEmpty()) {
            return;
        }

        int total = alignmentCounts.getTotalCount(pos);
        if (total == 0) {
            return;
        }

        Graphics2D graphics = context.getGraphics();
        for (String type : maCounts.getTypesInDrawOrder()) {
            int count = maCounts.getCount(type, pos);
            if (count == 0) {
                continue;
            }
            int height = Math.round((((float) count) / total) * barHeight);
            if (height <= 0) {
                continue;
            }
            graphics.setColor(MaColors.getColor(type));
            graphics.fillRect(pX, pBottom - height, dX, height);
        }
    }
}
