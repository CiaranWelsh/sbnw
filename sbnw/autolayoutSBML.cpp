/*== SAGITTARIUS =====================================================================
 * Copyright (c) 2012, Jesse K Medley
 * All rights reserved.

 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of The University of Washington nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//== BEGINNING OF CODE ===============================================================

//== INCLUDES ========================================================================

#include "autolayoutSBML.h"
#include "sbnw/error.h"

#include "sbml/SBMLTypes.h"
#include <sstream>

void gf_freeSBMLModel(gf_SBMLModel *lo = nullptr) {
    if (!lo)
        throw std::invalid_argument("Not a valid layout pointer");
    auto *doc = (libsbml::SBMLDocument *) lo->pdoc;
    delete doc;
    free(lo);
}

extern "C" gf_SBMLModel *gf_loadSBMLbuf(const char *buf) {
    auto *r = (gf_SBMLModel *) malloc(sizeof(gf_SBMLModel));
    libsbml::SBMLReader reader;
    libsbml::SBMLDocument *doc = reader.readSBMLFromString(buf);

    //not libSBML's documented way of failing, but just in case...
    if (!doc)
        throw std::invalid_argument("Failed to parse SBML");

    if (doc->getNumErrors()) {

        // if all are warnings, continue - else abort
        for (unsigned int i = 0; i < doc->getNumErrors(); ++i) {
            if (!doc->getError(i)->isWarning())
                return nullptr;
        }
    }

    r->pdoc = doc;
    return r;
}

extern "C" gf_SBMLModel *gf_loadSBMLfile(const char *path) {
    char *buf;
    size_t size = 0;
    FILE *file = nullptr;
    size_t bytes_read;


    file = fopen(path, "rb");
    if (!file)
        throw std::logic_error("Failed to open file, gf_loadSBMLfile");

    //get t3h s!z3
    fseek(file, 0, SEEK_END);
    size = ftell(file);
    rewind(file);

    assert(size > 0);

    //allocated buffer
    if (sizeof(char) != 1)
        throw std::logic_error("char must be one byte wide: gf_loadSBMLfile");
    buf = (char *) malloc(size + 1); //one extra byte for null char
    if (!buf)
        throw std::logic_error("Failed to allocate buffer: gf_loadSBMLfile");
    //read the whole file at once
    bytes_read = fread(buf, 1, size, file);
    if (bytes_read != size)
        throw std::logic_error("Failed to read whole file (wrong size specified?): gf_loadSBMLfile");
    //trip EOF indicator
    fgetc(file);
    if (!feof(file))
        throw std::logic_error("EOF Expected: gf_loadSBMLfile");
    buf[size] = '\0'; //terminating null char

    /*close*/
    fclose(file);

    gf_SBMLModel *mod = gf_loadSBMLbuf(buf);

    free(buf);

    return mod;

}



