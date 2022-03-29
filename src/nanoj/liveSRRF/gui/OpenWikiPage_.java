package nanoj.liveSRRF.gui;

import ij.plugin.PlugIn;
import static nanoj.liveSRRF.WebAuxiliary.OpenLink;

public class OpenWikiPage_ implements PlugIn {

    @Override
    public void run(String s) {

        String wiki_url = "https://github.com/HenriquesLab/NanoJ-eSRRF/wiki";
        OpenLink(wiki_url);

    }
}