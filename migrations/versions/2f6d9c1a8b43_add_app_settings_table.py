"""add app settings table

Revision ID: 2f6d9c1a8b43
Revises: b731f4c8d2a1
"""

from alembic import op
import sqlalchemy as sa


revision = "2f6d9c1a8b43"
down_revision = "b731f4c8d2a1"
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        "app_setting",
        sa.Column("key", sa.String(length=128), nullable=False),
        sa.Column("value", sa.Text(), nullable=False),
        sa.Column("description", sa.Text(), nullable=True),
        sa.Column("updated_at", sa.DateTime(), nullable=False),
        sa.PrimaryKeyConstraint("key"),
    )


def downgrade():
    op.drop_table("app_setting")
