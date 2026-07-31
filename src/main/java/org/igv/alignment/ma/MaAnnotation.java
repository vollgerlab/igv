package org.igv.alignment.ma;

/**
 * A single molecular annotation (MA tag) entry.
 * <p>
 * Coordinates are 0-based half-open and are given in molecular (original read) orientation,
 * i.e. counting from the leftmost base of the unaligned read.  The MA tag itself stores
 * 1-based closed positions; the conversion happens at parse time.
 *
 * @param type      annotation type name, e.g. "msp", "nuc", "fire"
 * @param strand    strand of this annotation, one of '+', '-', or '.'
 * @param start     0-based start position on the molecule
 * @param length    length in bases
 * @param qualities quality values from the AQ tag for this annotation, empty if the type has no quality spec
 * @param name      label from the AN tag, or null if unnamed
 */
public record MaAnnotation(String type, char strand, int start, int length, byte[] qualities, String name) {

    /**
     * Exclusive end position, 0-based.
     */
    public int end() {
        return start + length;
    }
}
