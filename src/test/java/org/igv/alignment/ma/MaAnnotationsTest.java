package org.igv.alignment.ma;

import org.junit.Test;

import java.util.List;

import static org.junit.Assert.*;

/**
 * Tests for the MA / AQ / AN tag parser.  The tag values used here are taken verbatim from the
 * molecular annotation spec examples.
 */
public class MaAnnotationsTest {

    @Test
    public void testNoQuality() {

        // MA:Z:1000;msp+:100-50,200-60
        MaAnnotations ma = MaAnnotations.parse("1000;msp+:100-50,200-60", null, null);

        assertNotNull(ma);
        assertEquals(1000, ma.readLength);
        assertEquals(List.of("msp"), ma.typeNames());
        assertEquals(2, ma.annotations.size());

        MaAnnotation first = ma.annotations.get(0);
        assertEquals("msp", first.type());
        assertEquals('+', first.strand());
        assertEquals(99, first.start());     // 1-based 100 -> 0-based 99
        assertEquals(50, first.length());
        assertEquals(149, first.end());
        assertEquals(0, first.qualities().length);
        assertNull(first.name());

        MaAnnotation second = ma.annotations.get(1);
        assertEquals(199, second.start());
        assertEquals(60, second.length());
        assertEquals(259, second.end());
    }

    @Test
    public void testPhredQuality() {

        // MA:Z:1000;msp+P:100-50,200-60  AQ:B:C,40,30
        byte[] aq = {40, 30};
        MaAnnotations ma = MaAnnotations.parse("1000;msp+P:100-50,200-60", aq, null);

        assertNotNull(ma);
        assertEquals(2, ma.annotations.size());
        assertEquals("P", ma.qualitySpec("msp"));
        assertArrayEquals(new byte[]{40}, ma.annotations.get(0).qualities());
        assertArrayEquals(new byte[]{30}, ma.annotations.get(1).qualities());
    }

    @Test
    public void testMultipleQualitiesPerAnnotation() {

        // MA:Z:1000;msp+PQ:100-50,200-60;nuc+:150-103   AQ:B:C,40,255,30,200
        byte[] aq = {40, (byte) 255, 30, (byte) 200};
        MaAnnotations ma = MaAnnotations.parse("1000;msp+PQ:100-50,200-60;nuc+:150-103", aq, null);

        assertNotNull(ma);
        assertEquals(List.of("msp", "nuc"), ma.typeNames());
        assertEquals(3, ma.annotations.size());

        // Quality values are grouped per annotation:  [msp1_P, msp1_Q, msp2_P, msp2_Q]
        assertArrayEquals(new byte[]{40, (byte) 255}, ma.annotations.get(0).qualities());
        assertArrayEquals(new byte[]{30, (byte) 200}, ma.annotations.get(1).qualities());

        // nuc has no quality spec, so it consumes nothing from AQ
        MaAnnotation nuc = ma.annotations.get(2);
        assertEquals("nuc", nuc.type());
        assertEquals(149, nuc.start());
        assertEquals(103, nuc.length());
        assertEquals(0, nuc.qualities().length);
        assertEquals("", ma.qualitySpec("nuc"));
    }

    @Test
    public void testPartialNames() {

        // MA:Z:1000;msp+P:100-50,200-60;nuc+:150-103,300-100   AQ:B:C,40,35   AN:Z:msp1,,,nuc2
        byte[] aq = {40, 35};
        MaAnnotations ma = MaAnnotations.parse(
                "1000;msp+P:100-50,200-60;nuc+:150-103,300-100", aq, "msp1,,,nuc2");

        assertNotNull(ma);
        assertEquals(4, ma.annotations.size());

        // AN is positional over all annotations, including the quality-less nuc entries
        assertEquals("msp1", ma.annotations.get(0).name());
        assertNull(ma.annotations.get(1).name());
        assertNull(ma.annotations.get(2).name());
        assertEquals("nuc2", ma.annotations.get(3).name());
    }

    @Test
    public void testRepeatedTypeNameAcrossStrands() {

        // MA:Z:10;ctcf+Q:1-4;ctcf-Q:6-3   AQ:B:C,200,180
        byte[] aq = {(byte) 200, (byte) 180};
        MaAnnotations ma = MaAnnotations.parse("10;ctcf+Q:1-4;ctcf-Q:6-3", aq, null);

        assertNotNull(ma);
        assertEquals(10, ma.readLength);
        assertEquals(List.of("ctcf"), ma.typeNames());
        assertEquals(2, ma.annotations.size());

        MaAnnotation forward = ma.annotations.get(0);
        assertEquals('+', forward.strand());
        assertEquals(0, forward.start());
        assertEquals(4, forward.end());
        assertArrayEquals(new byte[]{(byte) 200}, forward.qualities());

        MaAnnotation reverse = ma.annotations.get(1);
        assertEquals("ctcf", reverse.type());
        assertEquals('-', reverse.strand());
        assertEquals(5, reverse.start());
        assertEquals(8, reverse.end());
        assertArrayEquals(new byte[]{(byte) 180}, reverse.qualities());
    }

    @Test
    public void testToReadCoords() {

        MaAnnotations ma = MaAnnotations.parse("1000;msp+:100-50", null, null);
        assertNotNull(ma);
        MaAnnotation a = ma.annotations.get(0);   // start 99, length 50

        assertArrayEquals(new int[]{99, 149}, ma.toReadCoords(a, false));

        // Reverse strand alignments store the read reverse complemented:  [1000-149, 1000-99)
        assertArrayEquals(new int[]{851, 901}, ma.toReadCoords(a, true));
    }

    @Test
    public void testMalformed() {

        // A position without an inline length is not allowed
        assertNull(MaAnnotations.parse("1000;msp+:100", null, null));

        // Read length must be a number
        assertNull(MaAnnotations.parse("not-a-number;msp+:1-2", null, null));

        // AQ must hold enough values for every annotation of a quality bearing type
        assertNull(MaAnnotations.parse("1000;msp+P:100-50,200-60", new byte[]{40}, null));
        assertNull(MaAnnotations.parse("1000;msp+P:100-50", null, null));

        // Other malformations
        assertNull(MaAnnotations.parse(null, null, null));
        assertNull(MaAnnotations.parse("", null, null));
        assertNull(MaAnnotations.parse("1000;msp100-50", null, null));         // no strand indicator
        assertNull(MaAnnotations.parse("1000;msp+X:100-50", null, null));      // bad quality spec
        assertNull(MaAnnotations.parse("1000;msp+:0-50", null, null));         // positions are 1-based
        assertNull(MaAnnotations.parse("1000;msp+P:100-50;msp+Q:200-60",
                new byte[]{40, 30}, null));                                    // conflicting quality specs
    }
}
