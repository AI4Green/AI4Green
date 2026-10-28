import { formatDate, formatQuantity, itemState } from "./utils.jsx";

export function Badge({ tone, children }) {
  return <span className={`inv-badge inv-${tone}`}>{children}</span>;
}

const ACTIONS = [
  ["adjust", "Adjust"],
  ["move", "Move"],
  ["dispose", "Dispose"],
  ["history", "History"],
];

/**
 * One row per InventoryItem.
 * Pass showCompound for the flat "All containers" view, where items carry
 * compoundName and cas. Inside an expanded compound row it is omitted.
 */
export default function ContainersTable({
  items,
  showCompound = false,
  onItemAction,
}) {
  if (items.length === 0) {
    return <p className="inv-empty">No containers to show.</p>;
  }

  return (
    <table className="inv-table">
      <thead>
        <tr>
          {showCompound && <th>Compound</th>}
          <th>Location</th>
          <th>Quantity</th>
          <th>Lot</th>
          <th>Expiry</th>
          <th>State</th>
          <th>
            <span className="inv-sr">Actions</span>
          </th>
        </tr>
      </thead>
      <tbody>
        {items.map((item) => {
          const state = item.isArchived
            ? { label: "Archived", tone: "muted" }
            : itemState(item);
          const actions = item.isArchived
            ? ACTIONS.filter(([key]) => key === "history")
            : ACTIONS;

          return (
            <tr key={item.id}>
              {showCompound && (
                <td>
                  <b>{item.compoundName}</b>
                  <div className="inv-cas">{item.cas}</div>
                </td>
              )}
              <td>{item.location.name}</td>
              <td>
                {formatQuantity(item.quantity)} {item.unit}
              </td>
              <td>{item.lotNumber || "-"}</td>
              <td>{formatDate(item.expiryDate)}</td>
              <td>
                <Badge tone={state.tone}>{state.label}</Badge>
              </td>
              <td className="inv-actions">
                {actions.map(([key, label]) => (
                  <button
                    key={key}
                    type="button"
                    className="inv-link"
                    aria-label={`${label} container ${
                      item.lotNumber || item.id
                    }`}
                    onClick={() => onItemAction(key, item)}
                  >
                    {label}
                  </button>
                ))}
              </td>
            </tr>
          );
        })}
      </tbody>
    </table>
  );
}
