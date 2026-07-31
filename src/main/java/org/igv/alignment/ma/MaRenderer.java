package org.igv.alignment.ma;

import org.igv.alignment.*;

import java.awt.*;
import java.util.List;

/**
 * Draws MA (molecular annotation) intervals on an alignment row.
 * <p>
 * The analog of {@link org.igv.alignment.mods.BaseModificationRenderer}, except that an MA entry is
 * an interval (start + length) on the molecule rather than a single base.  Annotation coordinates
 * are lifted from read space to genome space through the alignment blocks, so an interval spanning
 * a deletion or an intron is drawn as one piece per block.
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

        // ponytail: this is a linear scan of every annotation on the read, on every repaint.  A
        // fiber-seq read carries thousands of nuc/msp entries, so the ceiling is
        // O(visible reads x annotations per read) per paint.  Annotations are not guaranteed to be
        // in ascending order, so there is no early exit on the outer loop; sorting them here would
        // cost more than the scan it saves.  The fix, if this ever shows up in a profile, is to sort
        // once at parse time in MaAnnotations and binary search the visible read-offset window.
        for (MaAnnotation a : annotations) {

            // Molecular coordinates -> offsets into the read sequence as stored in the BAM,
            // flipping for negative strand reads.
            int[] readCoords = ma.toReadCoords(a, negativeStrand);
            if (readCoords == null) {
                continue;
            }
            int readStart = readCoords[0];
            int readEnd = readCoords[1];
            if (readEnd <= readStart) {
                continue;
            }

            Color color = null;

            // Blocks are built in CIGAR order, so they ascend in both read offset and genome
            // position.  That allows the two breaks below.
            for (AlignmentBlock block : blocks) {

                // ponytail: reads stored without a sequence (SEQ = "*") have empty base arrays, so
                // there is nothing to map offsets against and their annotations are not drawn.
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

                // Intersect the annotation with the block in read offset space.  This is the CIGAR
                // liftover:  offset o within a block maps to genome position
                // block.getStart() + (o - blockStart).
                int fromOffset = Math.max(readStart, blockStart);
                int toOffset = Math.min(readEnd, blockEnd);

                int genomeStart = block.getStart() + (fromOffset - blockStart);
                int genomeEnd = block.getStart() + (toOffset - blockStart);

                int pX = (int) ((genomeStart - bpStart) / locScale);
                int pEnd = (int) ((genomeEnd - bpStart) / locScale);

                // Minimum width of 1 so short annotations survive zoom out.  Unlike the base
                // modification renderer these are not widened to 3 pixels -- they are intervals with
                // real boundaries, and padding them would misreport where a feature starts and ends.
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
