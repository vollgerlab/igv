package org.igv.alignment.ma;

import org.igv.alignment.AlignmentCounts;
import org.igv.alignment.RenderOptions;
import org.igv.track.RenderContext;

import java.awt.*;


/**
 * Draws MA (molecular annotation) composition into the coverage bar, the analog of
 * {@link org.igv.alignment.mods.BaseModificationCoverageRenderer}.
 * <p>
 * Types are overlaid, not stacked, each rising from the bottom at its own fraction of the covering
 * reads.  MA types are often nested rather than disjoint -- a FIRE is an MSP that cleared an FDR
 * threshold -- so stacking would report more than 100% of the reads at a position.
 */
public class MaCoverageRenderer {

    public static void drawAnnotations(RenderContext context,
                                       int pX,
                                       int pBottom,
                                       int dX,
                                       int barHeight,
                                       int pos,
                                       AlignmentCounts alignmentCounts,
                                       RenderOptions renderOptions) {

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
            if (renderOptions != null && !renderOptions.isMaTypeVisible(type)) {
                continue;
            }
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
