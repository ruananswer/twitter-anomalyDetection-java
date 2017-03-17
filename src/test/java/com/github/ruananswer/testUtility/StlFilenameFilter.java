package com.github.ruananswer.testUtility;

import java.io.File;
import java.io.FilenameFilter;

/**
 * Created by ruan on 16-4-18.
 */
public class StlFilenameFilter implements FilenameFilter {
    private String type;

    public StlFilenameFilter(String regex) {
        this.type = regex;
    }
    @Override
    public boolean accept(File dir, String name) {
        return name.endsWith(type);
    }
}
