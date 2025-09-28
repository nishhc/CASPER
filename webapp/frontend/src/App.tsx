import { BrowserRouter as Router, Routes, Route } from "react-router-dom";
import Home from "./Home";
import AssayDesignerPage from "./AssayDesignerPage";

export default function App() {
  return (
    <Router>
      <Routes>
        <Route path="/" element={<Home />} />
        <Route path="/assay" element={<AssayDesignerPage />} />
      </Routes>
    </Router>
  );
}
