with import <nixpkgs> { };
stdenv.mkDerivation {
  name = "pew";
  buildInputs = [ rustChannels.nightly.rust ];
}
