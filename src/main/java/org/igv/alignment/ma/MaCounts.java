package org.igv.alignment.ma;

import org.igv.alignment.Alignment;
import org.igv.alignment.AlignmentBlock;
import org.igv.alignment.ByteSubarray;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/**
 * Per position counts of MA (molecular annotation) intervals over a loaded alignment interval, the
 * analog of {@link org.igv.alignment.mods.BaseModificationCounts}.
 * <p>
 * Counts are dense arrays rather than the position maps the base modification counts use.  MA
 * intervals tile a molecule -- nucleosomes and MSPs together cover nearly every base of a fiber --
 * so a map would hold an entry per covered position anyway, at much greater cost per entry.
 */
public class MaCounts {

    private final int start;
    private final int end;

    /**
     * Annotation type name -> count of reads carrying that annotation at each position of
     * [start, end).  Arrays are allocated on first sighting of a type, so a file with no MA tags
     * costs nothing beyond the empty map.  Insertion ordered so the coverage renderer sees types in
     * the order the data presents them.
     */
    private final Map<String, int[]> countsByType = new LinkedHashMap<>();

    /**
     * Types ordered back to front, see {@link #getTypesInDrawOrder()}.  Rebuilt whenever a new type
     * turns up, which in practice means once per interval.
     */
    private List<String> drawOrder;

    public MaCounts(int start, int end) {
        this.start = start;
        this.end = end;
    }

    /**
     * Increment counts for every genome position covered by an MA annotation on this alignment.
     *
     * @param alignment alignment to tally
     */
    public void incrementCounts(Alignment alignment) {

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

        for (MaAnnotation a : annotations) {

            int[] readCoords = ma.toReadCoords(a, negativeStrand);
            if (readCoords == null) {
                continue;
            }
            int readStart = readCoords[0];
            int readEnd = readCoords[1];
            if (readEnd <= readStart) {
                continue;
            }

            int[] counts = null;

            // Same liftover MaRenderer performs, tallying positions instead of filling pixels.
            // Blocks ascend in both read offset and genome position, hence the breaks.
            for (AlignmentBlock block : blocks) {

                if (block.isSoftClip()) {
                    continue;           // matches BaseModificationCounts, which ignores soft clips
                }

                ByteSubarray bases = block.getBases();
                if (bases == null || bases.length == 0) {
                    continue;
                }

                int blockStart = bases.startOffset;
                int blockEnd = blockStart + bases.length;

                if (blockStart >= readEnd) {
                    break;
                }
                if (blockEnd <= readStart) {
                    continue;
                }

                int fromOffset = Math.max(readStart, blockStart);
                int toOffset = Math.min(readEnd, blockEnd);

                int genomeStart = block.getStart() + (fromOffset - blockStart);
                int genomeEnd = block.getStart() + (toOffset - blockStart);

                // Clip to the interval this counts object covers
                int from = Math.max(genomeStart, start);
                int to = Math.min(genomeEnd, end);
                if (to <= from) {
                    continue;
                }

                if (counts == null) {
                    counts = countsByType.computeIfAbsent(a.type(), t -> new int[end - start]);
                }
                for (int p = from; p < to; p++) {
                    counts[p - start]++;
                }
            }
        }
    }

    /**
     * @return annotation type names seen in this interval, in first seen order
     */
    public Set<String> getTypes() {
        return Collections.unmodifiableSet(countsByType.keySet());
    }

    /**
     * Types ordered back to front for drawing:  nucleosomes behind, then anything unrecognized, then
     * MSPs, with FIREs painted last so they always read along the bottom of the coverage bar.
     * <p>
     * The order is fixed rather than derived from the counts at each position.  Ordering by count
     * flips the z order from column to column wherever two types cross over, which makes a type
     * appear to jump between foreground and background across the track.
     */
    public List<String> getTypesInDrawOrder() {
        if (drawOrder == null || drawOrder.size() != countsByType.size()) {
            List<String> types = new ArrayList<>(countsByType.keySet());
            // Stable, so unrecognized types keep first seen order among themselves
            types.sort(Comparator.comparingInt(MaCounts::drawRank));
            drawOrder = types;
        }
        return drawOrder;
    }

    private static int drawRank(String type) {
        switch (type.toLowerCase(Locale.ROOT)) {
            case "nuc":
                return 0;
            case "msp":
                return 2;
            case "fire":
                return 3;
            default:
                return 1;
        }
    }

    /**
     * @param type annotation type name
     * @param pos  genome position
     * @return number of reads carrying an annotation of this type over this position
     */
    public int getCount(String type, int pos) {
        int[] counts = countsByType.get(type);
        if (counts == null || pos < start || pos >= end) {
            return 0;
        }
        return counts[pos - start];
    }

    public boolean isEmpty() {
        return countsByType.isEmpty();
    }
}
