package org.igv.alignment.ma;

/**
 * A single molecular annotation (MA tag) entry.  Coordinates are 0-based half-open in molecular
 * (original read) orientation; the tag stores 1-based closed positions and is converted at parse
 * time.  Strand is '+', '-', or '.'.
 */
public record MaAnnotation(String type, char strand, int start, int length, byte[] qualities, String name) {

    public int end() {
        return start + length;
    }
}
