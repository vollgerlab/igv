package org.igv.alignment.ma;

import org.igv.alignment.*;

import java.awt.*;
import java.util.List;

/**
 * Draws MA (molecular annotation) intervals on an alignment row.  The analog of
 * {@link org.igv.alignment.mods.BaseModificationRenderer}, except an MA entry is an interval rather
 * than a single base.  Coordinates are lifted through the alignment blocks, so an interval spanning
 * a deletion or intron is drawn as one piece per block.
 */
public class MaRenderer {

    public static void drawAnnotations(
            Alignment alignment,
            double bpStart,
            double locScale,
            Rectangle rowRect,
            Graphics g,
            RenderOptions renderOptions) {

        MaAnnotations ma = MaAnnotations.forAlignment(alignment);
        if (ma == null) {
            return;
        }

        List<MaAnnotation> annotations = ma.annotations;
        if (annotations == null || annotations.isEmpty()) {
            return;
        }

        AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
        if (blocks == null || blocks.length == 0) {
            return;
        }

        boolean negativeStrand = alignment.isNegativeStrand();

        // Determine if we should leave a margin (same logic as AlignmentRenderer)
        boolean leaveMargin = rowRect.height > 2;

        int pY = (int) rowRect.getY();
        int rectHeight = Math.max(1, rowRect.height - (leaveMargin ? 2 : 0));
        int clipLeft = (int) rowRect.getX();
        int clipRight = (int) rowRect.getMaxX();

        // Linear scan per repaint, O(visible reads x annotations per read).  Annotations are not
        // guaranteed sorted, so there is no early exit here; sorting would cost more than it saves.
        for (MaAnnotation a : annotations) {

            // Molecular coordinates -> offsets into the read sequence as stored in the BAM
            int[] readCoords = ma.toReadCoords(a, negativeStrand);
            if (readCoords == null) {
                continue;
            }
            int readStart = readCoords[0];
            int readEnd = readCoords[1];
            if (readEnd <= readStart) {
                continue;
            }

            if (renderOptions != null && !renderOptions.isMaTypeVisible(a.type())) {
                continue;
            }

            Color color = null;

            // Blocks ascend in read offset and genome position, hence the breaks below
            for (AlignmentBlock block : blocks) {

                // Reads stored without a sequence (SEQ = "*") have nothing to map offsets against
                ByteSubarray bases = block.getBases();
                if (bases == null || bases.length == 0) {
                    continue;
                }

                int blockStart = bases.startOffset;
                int blockEnd = blockStart + bases.length;

                if (blockStart >= readEnd) {
                    break;              // past the annotation, and blocks only go right from here
                }
                if (blockEnd <= readStart) {
                    continue;           // not there yet
                }

                // CIGAR liftover: offset o in a block maps to block.getStart() + (o - blockStart)
                int fromOffset = Math.max(readStart, blockStart);
                int toOffset = Math.min(readEnd, blockEnd);

                int genomeStart = block.getStart() + (fromOffset - blockStart);
                int genomeEnd = block.getStart() + (toOffset - blockStart);

                int pX = (int) ((genomeStart - bpStart) / locScale);
                int pEnd = (int) ((genomeEnd - bpStart) / locScale);

                // Minimum width 1 so short annotations survive zoom out.  Not widened to 3px like
                // the basemod renderer -- these have real boundaries and padding would misreport them.
                int dX = Math.max(1, pEnd - pX);

                // Don't draw out of the clipping rect
                if (pX > clipRight) {
                    break;
                } else if (pX + dX < clipLeft) {
                    continue;
                }

                if (color == null) {
                    color = MaColors.getColor(a);
                }
                g.setColor(color);
                g.fillRect(pX, pY, dX, rectHeight);
            }
        }
    }
}
