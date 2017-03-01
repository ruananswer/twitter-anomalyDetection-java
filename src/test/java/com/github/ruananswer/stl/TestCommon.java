package com.github.ruananswer.stl;

import org.apache.commons.io.IOUtils;

import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;

/**
 * Created by ruan on 16-4-18.
 */
public class TestCommon {
   public static String getResourceAsString(String resource) throws IOException {
        InputStream is = ClassLoader.getSystemResourceAsStream(resource);
        StringWriter writer = new StringWriter();
        IOUtils.copy(is, writer);
        return writer.toString();
    }
}
