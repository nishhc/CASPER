import React from "react";
import "./home.css";

export default function Home() {
  return (
    <section id="hero" className="hero">
      {/* background contour lines */}


      <div className="content">
        <h1 className="title">CASPER</h1>
        <p className="subtitle">Integrated RPA + CRISPR-Cas12 Assay Designer</p>
        <div className="buttons">
          <a href="#about" className="pill">About</a>
          <a href="/assay" className="pill">Generate</a>
          <a href="#rank" className="pill">Rank Custom</a>
        </div>
      </div>
    </section>
  );
}
