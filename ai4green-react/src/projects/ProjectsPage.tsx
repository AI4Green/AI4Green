import React from "react";

export interface Props {
  name: string;
  enthusiasmLevel?: number;
}

function Hello({ name, enthusiasmLevel = 1 }: Props) {
  if (enthusiasmLevel <= 0) {
    throw new Error("You should be more enthusiastic");
  }

  return (
    <div className="hello">
      <div className="greeting">
        Hello {name + getExclamationMarks(enthusiasmLevel)}
      </div>
    </div>
  );
}

// helpers

function getExclamationMarks(numChars: number) {
  return Array(numChars + 1).join("!!!");
}

class helloClass extends React.component<Props> {
  render() {
    const { name, enthusiasmLevel = 1 } = this.Props;

    if (enthusiasmLevel <= 0) {
      throw new Error("Be more enthusiastic");
    }

    return (
      <div className="hello">
        <div className="greeting">
          Hello {name + getExclamationMarks(enthusiasmLevel)}
        </div>
      </div>
    );
  }
}

export default helloClass;
