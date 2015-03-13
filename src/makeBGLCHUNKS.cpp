#include "utils.h"

int main(int argc, char **argv) {
  string buffer;

  char *vcf = NULL;
  int window = -1;
  int overlap = -1;
  char *output = NULL;
  bool hardChunkSizeLimit = false;
  bool chunksEqualSize = false;
  bool noOutOffset = false;

  // Reading arguments
  for (int a = 1; a < argc; a++) {
    if (strcmp(argv[a], "--vcf") == 0) {
      vcf = argv[a + 1];
    }
    if (strcmp(argv[a], "--window") == 0) {
      window = atoi(argv[a + 1]);
    }
    if (strcmp(argv[a], "--overlap") == 0) {
      overlap = atoi(argv[a + 1]);
    }
    if (strcmp(argv[a], "--output") == 0) {
      output = argv[a + 1];
    }
    if (strcmp(argv[a], "--hard_chunk_size_limit") == 0) {
      hardChunkSizeLimit = true;
      cout << endl
           << "Creating chunks with an upper limit of 2*overlap + window sites."
           << endl;
    }
    if (strcmp(argv[a], "--make_chunks_equal_size") == 0) {
      chunksEqualSize = true;
      hardChunkSizeLimit = true;
      cout << endl << "Creating chunks with exactly 2*overlap + window sites.  "
                      "Second chunk will start after window + overlap sites."
           << endl;
    }
    if (strcmp(argv[a], "--make_chunks_equal_size") == 0) {
      chunksEqualSize = true;
      hardChunkSizeLimit = true;
      cout << endl << "Creating chunks with exactly 2*overlap + window sites.  "
                      "Second chunk will start after window + overlap sites."
           << endl;
    }
    if (strcmp(argv[a], "--no_out_offset") == 0) {
      noOutOffset = true;
      cout << endl << "Output chunk sites are inclusive" << endl;
    }
  }

  // Checking arguments
  if (vcf == NULL) {
    cout << "Specify a VCF file with --vcf" << endl;
    exit(0);
  } else
    cout << "VCF:\t" << vcf << endl;
  if (window < 0) {
    cout << "Specify a window size in number of sites with --window" << endl;
    exit(0);
  } else
    cout << "Window:\t" << window << endl;
  if (overlap < 0) {
    cout << "Specify an overlap size in number of sites with --overlap" << endl;
    exit(0);
  } else
    cout << "Overlap:\t" << overlap << endl;
  if (output == NULL) {
    cout << "Specify an output file with --output" << endl;
    exit(0);
  } else
    cout << "Output:\t" << output << endl;

  // Reading VCF file
  string chr;
  vector<int> P;
  cout << "Reading sites in [" << vcf << "]" << endl;
  ifile fds(vcf);
  size_t lineNum = 0;
  vector<string> tokens;
  while (getline(fds, buffer)) {
    if (buffer[0] != '#') {
      tokens = sutils::tokenize(buffer, "\t", 5);
      if (chr.empty())
        chr = tokens[0];
      P.push_back(atoi(tokens[1].c_str()));
      ++lineNum;
    }
  }
  fds.close();
  cout << "  * " << P.size() << " sites read" << endl;
  double density = (P.back() - P[0]) * 1.0 / P.size();
  cout << "  * " << sutils::double2str(density, 1) << " bp betweem sites"
       << endl;
  cout << "  * " << sutils::double2str(density * window, 1) << " bp per chunk"
       << endl;
  cout << "  * " << sutils::double2str(density * overlap, 1)
       << " bp per overlap" << endl;

  // Building the chunks
  cout << "Building chunks" << endl;
  vector<int> B, E;
  for (int s = 0; s < P.size(); s += window) {
    B.push_back(s);
    if (s == 0 && chunksEqualSize)
      s += overlap;
    E.push_back(s + window - 1);
  }
  if (E.back() >= P.size())
    E.back() = P.size() - 1;

  if (hardChunkSizeLimit) {
    B.back() = P.size() - window - overlap;
    E.back() = P.size() - overlap - 1;
  } else if ((E.back() - B.back()) < (window / 2)) {
    B.erase(B.begin() + B.size() - 1);
    E.erase(E.begin() + E.size() - 1);
    E.back() = P.size() - 1;
  }

  for (int w = 0; w < B.size(); w++) {
    if (B[w] > 0)
      B[w] -= overlap;
    if (E[w] < P.size() - 1)
      E[w] += overlap;
  }

  for (int w = 0; w < B.size(); w++) {
    if (noOutOffset) {
        B[w] = P.size() > B[w] ? P.at(B[w]) : P.back();
        E[w] = P.size() > E[w] ? P.at(E[w]) : P.back();
    } else {
      B[w] = P[B[w]] - 1;
      E[w] = P[E[w]] + 1;
    }
  }
  cout << "  * " << B.size() << " chunks built" << endl;

  cerr << "Writing chunks in [" << output << "]" << endl;
  ofile fdo(output);
  for (int w = 0; w < B.size(); w++) {
    fdo << chr << "\t" << B[w] << "\t" << E[w] << endl;
  }
  fdo.close();
}
