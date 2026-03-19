import React from "react";
import "./App.css";
import helloClass from "./projects/ProjectsPage";

function App() {
  return (
    <div className="App">
      <helloClass name="David" enthusiasmLevel={5}></helloClass>
    </div>
  );
}

export default App;
