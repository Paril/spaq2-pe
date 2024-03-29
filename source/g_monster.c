/*
Copyright (C) 1997-2001 Id Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/
#include "g_local.h"


//
// monster weapons
//

//FIXME mosnters should call these with a totally accurate direction
// and we can mess it up based on skill.  Spread should be for normal
// and we can tighten or loosen based on skill.  We could muck with
// the damages too, but I'm not sure that's such a good idea.
// SPAQ
// doing exactly this!!
static void monster_setup_accuracy(edict_t *self, vec3_t start, float **dir)
{
    edict_t *enemy = self->enemy;

    // ???
    if (!enemy || *dir)
        return;

    static vec3_t vdir;
    *dir = &vdir[0];

    if (gi.trace(start, NULL, NULL, start, self, MASK_SHOT).fraction < 1.0)
    {
        AngleVectors(self->s.angles, vdir, NULL, NULL);
        return;
    }
    
    vec3_t point;
    float dist;
    VectorSubtract(self->s.origin, self->enemy->s.origin, point);
    dist = VectorLength(point);

    // pick a random point on the enemy's bbox to aim at.
    VectorAdd(self->enemy->absmin, self->enemy->absmax, point);
    VectorScale(point, 0.5f, point);

    for (int i = 0; i < 3; i++)
        point[i] += (crandom() * (self->enemy->size[i] * 0.48f + (crandom() * (3.0 - skill->value) * self->enemy->size[i])));

    // 0-1 for hardness (0 is easy, 1 is hard)
    float hardness = skill->value / 3.f;

    // miss; factor in velocity.
    // the miss factor varies upon velocity
    float r = crandom() * max(0.25f, 1.0 - hardness);
    VectorMA(point, FRAMETIME * r, self->enemy->velocity, point);

    VectorSubtract(point, start, vdir);
}
// SPAQ
void monster_fire_bullet(edict_t *self, vec3_t start, vec3_t dir, int damage, int kick, int hspread, int vspread, int flashtype)
{
    monster_setup_accuracy(self, start, (float **) &dir);

    if (MONSTER_IS_CHAMP(self, CHAMPION_STRENGTH))
    {
        fire_bullet_sparks(self, start, dir, damage * 2, kick, 0, 0, MOD_M4 | MOD_MONSTER);
	    gi.sound(self, CHAN_WEAPON, gi.soundindex("weapons/m4a1fire.wav"), 1, ATTN_LOUD, 0);

		EjectShell(self, start, MOD_M4);
    }
    else
    {
        fire_bullet(self, start, dir, damage, kick, 0, 0, MOD_MP5 | MOD_MONSTER);
	    gi.sound(self, CHAN_WEAPON, gi.soundindex("weapons/mp5fire1.wav"), 1, ATTN_LOUD, 0);

		EjectShell(self, start, MOD_MP5);
    }
    
    gi.WriteByte(svc_muzzleflash2);
    gi.WriteShort(self - g_edicts);
    gi.WriteByte(flashtype);
    gi.multicast(start, MULTICAST_PVS);
}

void monster_fire_shotgun(edict_t *self, vec3_t start, vec3_t aimdir, int damage, int kick, int hspread, int vspread, int count, int flashtype)
{
    monster_setup_accuracy(self, start, (float **) &aimdir);

    if (MONSTER_IS_CHAMP(self, CHAMPION_STRENGTH))
    {
        fire_shotgun(self, start, aimdir, damage * 2, kick * 2, hspread * 2, vspread * 2, count * 2, MOD_HC | MOD_MONSTER);
	    gi.sound(self, CHAN_WEAPON, gi.soundindex("weapons/cannon_fire.wav"), 1, ATTN_LOUD, 0);
		
		EjectShell(self, start, MOD_HC);
		EjectShell(self, start, MOD_HC);
    }
    else
    {
        fire_shotgun(self, start, aimdir, damage, kick, hspread, vspread, count, MOD_M3 | MOD_MONSTER);
	    gi.sound(self, CHAN_WEAPON, gi.soundindex("weapons/shotgf1b.wav"), 1, ATTN_LOUD, 0);

		EjectShell(self, start, MOD_M3);
    }

    gi.WriteByte(svc_muzzleflash2);
    gi.WriteShort(self - g_edicts);
    gi.WriteByte(flashtype);
    gi.multicast(start, MULTICAST_PVS);
}

void monster_fire_blaster(edict_t *self, vec3_t start, vec3_t dir, int damage, int speed, int flashtype, int effect)
{
    //fire_blaster(self, start, dir, damage, speed, effect, false);

    if (damage >= 30)
    {
        monster_setup_accuracy(self, start, (float **) &dir);
        fire_bullet_sniper(self, start, dir, damage, damage, 0, 0, MOD_SNIPER | MOD_MONSTER);
        gi.sound(self, CHAN_WEAPON, gi.soundindex("weapons/ssgfire.wav"), 1, ATTN_LOUD, 0);

		EjectShell(self, start, MOD_SNIPER);
    }
    else
    {
        for (int i = 0; i < (MONSTER_IS_CHAMP(self, CHAMPION_STRENGTH) ? 2 : 1); i++)
        {
            float *r = dir;
            monster_setup_accuracy(self, start, (float **) &dir);
            fire_bullet(self, start, dir, damage, damage, 0, 0, MOD_MK23 | MOD_MONSTER);
            gi.sound(self, CHAN_WEAPON, gi.soundindex("weapons/mk23fire.wav"), 1, ATTN_LOUD, (float) i * 0.06f);
            dir = r;
			EjectShell(self, start, MOD_MK23);
        }
    }

    gi.WriteByte(svc_muzzleflash2);
    gi.WriteShort(self - g_edicts);
    gi.WriteByte(flashtype);
    gi.multicast(start, MULTICAST_PVS);
}

void monster_fire_grenade(edict_t *self, vec3_t start, vec3_t aimdir, int damage, int speed, int flashtype)
{
    monster_setup_accuracy(self, start, (float **) &aimdir);

    fire_grenade2(self, start, aimdir, damage * 3, speed, 2 * HZ, damage * 6, false);

    if (MONSTER_IS_CHAMP(self, CHAMPION_STRENGTH))
        fire_grenade2(self, start, aimdir, damage, speed * 2, 0.8 * HZ, damage * 4, false);

    gi.WriteByte(svc_muzzleflash2);
    gi.WriteShort(self - g_edicts);
    gi.WriteByte(flashtype);
    gi.multicast(start, MULTICAST_PVS);
}

void monster_fire_rocket(edict_t *self, vec3_t start, vec3_t dir, int damage, int speed, int flashtype)
{
    monster_setup_accuracy(self, start, (float **) &dir);

    fire_grenade2(self, start, dir, damage * 4, speed, 2 * HZ, damage * 8, false);

    if (MONSTER_IS_CHAMP(self, CHAMPION_STRENGTH))
        fire_grenade2(self, start, dir, damage * 2, speed * 2, 0.8 * HZ, damage * 6, false);

    gi.WriteByte(svc_muzzleflash2);
    gi.WriteShort(self - g_edicts);
    gi.WriteByte(flashtype);
    gi.multicast(start, MULTICAST_PVS);
}

void monster_fire_railgun(edict_t *self, vec3_t start, vec3_t aimdir, int damage, int kick, int flashtype)
{
    monster_setup_accuracy(self, start, (float **) &aimdir);

    fire_bullet_sniper(self, start, aimdir, damage, damage, 0, 0, MOD_SNIPER | MOD_MONSTER);
    gi.sound(self, CHAN_WEAPON, gi.soundindex("weapons/ssgfire.wav"), 1, ATTN_LOUD, 0);
	EjectShell(self, start, MOD_SNIPER);

    gi.WriteByte(svc_muzzleflash2);
    gi.WriteShort(self - g_edicts);
    gi.WriteByte(flashtype);
    gi.multicast(start, MULTICAST_PVS);
}

void monster_fire_bfg(edict_t *self, vec3_t start, vec3_t aimdir, int damage, int speed, int kick, float damage_radius, int flashtype)
{
    monster_setup_accuracy(self, start, (float **) &aimdir);

    fire_bullet_sniper(self, start, aimdir, damage, damage, 0, 0, MOD_SNIPER | MOD_MONSTER);
    gi.sound(self, CHAN_WEAPON, gi.soundindex("weapons/ssgfire.wav"), 1, ATTN_LOUD, 0);
	EjectShell(self, start, MOD_SNIPER);

    gi.WriteByte(svc_muzzleflash2);
    gi.WriteShort(self - g_edicts);
    gi.WriteByte(flashtype);
    gi.multicast(start, MULTICAST_PVS);
}



//
// Monster utility functions
//

void M_FliesOff(edict_t *self)
{
    self->s.effects &= ~EF_FLIES;
    self->s.sound = 0;
}

void M_FliesOn(edict_t *self)
{
    if (self->waterlevel)
        return;
    self->s.effects |= EF_FLIES;
    self->s.sound = gi.soundindex("infantry/inflies1.wav");
    self->think = M_FliesOff;
    self->nextthink = level.framenum + 60 * HZ;
}

void M_FlyCheck(edict_t *self)
{
    if (self->waterlevel)
        return;

    if (random() > 0.5f)
        return;

    self->think = M_FliesOn;
    self->nextthink = level.framenum + (5 + 10 * random()) * HZ;
}

void AttackFinished(edict_t *self, float time)
{
    self->monsterinfo.attack_finished = level.framenum + time * HZ;
}


void M_CheckGround(edict_t *ent)
{
    vec3_t      point;
    trace_t     trace;

    if (ent->flags & (FL_SWIM | FL_FLY))
        return;

    if (ent->velocity[2] > 100) {
        ent->groundentity = NULL;
        return;
    }

// if the hull point one-quarter unit down is solid the entity is on ground
    point[0] = ent->s.origin[0];
    point[1] = ent->s.origin[1];
    point[2] = ent->s.origin[2] - 0.25f;

    trace = gi.trace(ent->s.origin, ent->mins, ent->maxs, point, ent, MASK_MONSTERSOLID);

    // check steepness
    if (trace.plane.normal[2] < 0.7f && !trace.startsolid) {
        ent->groundentity = NULL;
        return;
    }

//  ent->groundentity = trace.ent;
//  ent->groundentity_linkcount = trace.ent->linkcount;
//  if (!trace.startsolid && !trace.allsolid)
//      VectorCopy (trace.endpos, ent->s.origin);
    if (!trace.startsolid && !trace.allsolid) {
        VectorCopy(trace.endpos, ent->s.origin);
        ent->groundentity = trace.ent;
        ent->groundentity_linkcount = trace.ent->linkcount;
        ent->velocity[2] = 0;
    }
}


void M_CatagorizePosition(edict_t *ent)
{
    vec3_t      point;
    int         cont;

//
// get waterlevel
//
    point[0] = ent->s.origin[0];
    point[1] = ent->s.origin[1];
    point[2] = ent->s.origin[2] + ent->mins[2] + 1;
    cont = gi.pointcontents(point);

    if (!(cont & MASK_WATER)) {
        ent->waterlevel = 0;
        ent->watertype = 0;
        return;
    }

    ent->watertype = cont;
    ent->waterlevel = 1;
    point[2] += 26;
    cont = gi.pointcontents(point);
    if (!(cont & MASK_WATER))
        return;

    ent->waterlevel = 2;
    point[2] += 22;
    cont = gi.pointcontents(point);
    if (cont & MASK_WATER)
        ent->waterlevel = 3;
}


void M_WorldEffects(edict_t *ent)
{
    int     dmg;

    if (ent->health > 0) {
        if (!(ent->flags & FL_SWIM)) {
            if (ent->waterlevel < 3) {
                ent->air_finished_framenum = level.framenum + 12 * HZ;
            } else if (ent->air_finished_framenum < level.framenum) {
                // drown!
                if (ent->pain_debounce_framenum < level.framenum) {
                    dmg = 2 + 2 * ((level.framenum - ent->air_finished_framenum) / HZ);
                    if (dmg > 15)
                        dmg = 15;
                    T_Damage(ent, world, world, vec3_origin, ent->s.origin, vec3_origin, dmg, 0, DAMAGE_NO_ARMOR, MOD_WATER);
                    ent->pain_debounce_framenum = level.framenum + 1 * HZ;
                }
            }
        } else {
            if (ent->waterlevel > 0) {
                ent->air_finished_framenum = level.framenum + 9 * HZ;
            } else if (ent->air_finished_framenum < level.framenum) {
                // suffocate!
                if (ent->pain_debounce_framenum < level.framenum) {
                    dmg = 2 + 2 * ((level.framenum - ent->air_finished_framenum) / HZ);
                    if (dmg > 15)
                        dmg = 15;
                    T_Damage(ent, world, world, vec3_origin, ent->s.origin, vec3_origin, dmg, 0, DAMAGE_NO_ARMOR, MOD_WATER);
                    ent->pain_debounce_framenum = level.framenum + 1 * HZ;
                }
            }
        }
    }

    if (ent->waterlevel == 0) {
        if (ent->flags & FL_INWATER) {
            gi.sound(ent, CHAN_BODY, gi.soundindex("player/watr_out.wav"), 1, ATTN_NORM, 0);
            ent->flags &= ~FL_INWATER;
        }
        return;
    }

    if ((ent->watertype & CONTENTS_LAVA) && !(ent->flags & FL_IMMUNE_LAVA)) {
        if (ent->damage_debounce_framenum < level.framenum) {
            ent->damage_debounce_framenum = level.framenum + 0.2f * HZ;
            T_Damage(ent, world, world, vec3_origin, ent->s.origin, vec3_origin, 10 * ent->waterlevel, 0, 0, MOD_LAVA);
        }
    }
    if ((ent->watertype & CONTENTS_SLIME) && !(ent->flags & FL_IMMUNE_SLIME)) {
        if (ent->damage_debounce_framenum < level.framenum) {
            ent->damage_debounce_framenum = level.framenum + 1 * HZ;
            T_Damage(ent, world, world, vec3_origin, ent->s.origin, vec3_origin, 4 * ent->waterlevel, 0, 0, MOD_SLIME);
        }
    }

    if (!(ent->flags & FL_INWATER)) {
        if (!(ent->svflags & SVF_DEADMONSTER)) {
            if (ent->watertype & CONTENTS_LAVA)
                if (random() <= 0.5f)
                    gi.sound(ent, CHAN_BODY, gi.soundindex("player/lava1.wav"), 1, ATTN_NORM, 0);
                else
                    gi.sound(ent, CHAN_BODY, gi.soundindex("player/lava2.wav"), 1, ATTN_NORM, 0);
            else if (ent->watertype & CONTENTS_SLIME)
                gi.sound(ent, CHAN_BODY, gi.soundindex("player/watr_in.wav"), 1, ATTN_NORM, 0);
            else if (ent->watertype & CONTENTS_WATER)
                gi.sound(ent, CHAN_BODY, gi.soundindex("player/watr_in.wav"), 1, ATTN_NORM, 0);
        }

        ent->flags |= FL_INWATER;
        ent->damage_debounce_framenum = 0;
    }
}


void M_droptofloor(edict_t *ent)
{
    vec3_t      end;
    trace_t     trace;

    ent->s.origin[2] += 1;
    VectorCopy(ent->s.origin, end);
    end[2] -= 256;

    trace = gi.trace(ent->s.origin, ent->mins, ent->maxs, end, ent, MASK_MONSTERSOLID);

    if (trace.fraction == 1 || trace.allsolid)
        return;

    VectorCopy(trace.endpos, ent->s.origin);

    gi.linkentity(ent);
    M_CheckGround(ent);
    M_CatagorizePosition(ent);
}


void M_SetEffects(edict_t *ent)
{
    ent->s.effects &= ~(EF_COLOR_SHELL | EF_POWERSCREEN);
    ent->s.renderfx &= ~(RF_SHELL_RED | RF_SHELL_GREEN | RF_SHELL_BLUE);

    if (ent->monsterinfo.aiflags & AI_RESURRECTING) {
        ent->s.effects |= EF_COLOR_SHELL;
        ent->s.renderfx |= RF_SHELL_RED;
    }

    if (ent->health <= 0)
        return;

    if (ent->powerarmor_framenum > level.framenum) {
        if (ent->monsterinfo.power_armor_type == POWER_ARMOR_SCREEN) {
            ent->s.effects |= EF_POWERSCREEN;
        } else if (ent->monsterinfo.power_armor_type == POWER_ARMOR_SHIELD) {
            ent->s.effects |= EF_COLOR_SHELL;
            ent->s.renderfx |= RF_SHELL_GREEN;
        }
    }

    // SPAQ
    if (MONSTER_IS_CHAMP(ent, CHAMPION_HASTE))
        ent->s.effects |= EF_BLASTER | EF_IONRIPPER;
    if (MONSTER_IS_CHAMP(ent, CHAMPION_SLOW_BLEED))
        ent->s.effects |= EF_GIB;
    if (MONSTER_IS_CHAMP(ent, CHAMPION_STRENGTH))
        ent->s.effects |= EF_QUAD;
    if (MONSTER_IS_CHAMP(ent, CHAMPION_STRONK))
        ent->s.effects |= EF_PENT;
    if (MONSTER_IS_CHAMP(ent, CHAMPION_VEST) || MONSTER_IS_CHAMP(ent, CHAMPION_HELMET))
        ent->s.effects |= EF_DOUBLE;
    // SPAQ
}


void M_MoveFrame(edict_t *self)
{
    mmove_t *move;
    int     index;

    move = self->monsterinfo.currentmove;
    self->nextthink = level.framenum + 1;

    if ((self->monsterinfo.nextframe) && (self->monsterinfo.nextframe >= move->firstframe) && (self->monsterinfo.nextframe <= move->lastframe)) {
        self->s.frame = self->monsterinfo.nextframe;
        self->monsterinfo.nextframe = 0;
    } else {
        if (self->s.frame == move->lastframe) {
            if (move->endfunc) {
                move->endfunc(self);

                // regrab move, endfunc is very likely to change it
                move = self->monsterinfo.currentmove;

                // check for death
                if (self->svflags & SVF_DEADMONSTER)
                    return;
            }
        }

        if (self->s.frame < move->firstframe || self->s.frame > move->lastframe) {
            self->monsterinfo.aiflags &= ~AI_HOLD_FRAME;
            self->s.frame = move->firstframe;
        } else {
            if (!(self->monsterinfo.aiflags & AI_HOLD_FRAME)) {
                self->s.frame++;
                if (self->s.frame > move->lastframe)
                    self->s.frame = move->firstframe;
            }
        }
    }

    index = self->s.frame - move->firstframe;
    if (move->frame[index].aifunc) {
        if (!(self->monsterinfo.aiflags & AI_HOLD_FRAME))
            move->frame[index].aifunc(self, move->frame[index].dist * self->monsterinfo.scale);
        else
            move->frame[index].aifunc(self, 0);
    }

    if (move->frame[index].thinkfunc)
        move->frame[index].thinkfunc(self);
}

// SPAQ
void Do_Bleeding (edict_t *ent);
// SPAQ

void monster_think(edict_t *self)
{
    M_MoveFrame(self);
    // SPAQ
    if ((self->monsterinfo.champion & CHAMPION_HASTE) && (level.framenum % 3) == 0)
        M_MoveFrame(self);
    // SPAQ
    if (self->linkcount != self->monsterinfo.linkcount) {
        self->monsterinfo.linkcount = self->linkcount;
        M_CheckGround(self);
    }
    M_CatagorizePosition(self);
    M_WorldEffects(self);
    M_SetEffects(self);

    // SPAQ
	// zucc handle any bleeding damage here
	Do_Bleeding(self);
    // SPAQ
}


/*
================
monster_use

Using a monster makes it angry at the current activator
================
*/
void monster_use(edict_t *self, edict_t *other, edict_t *activator)
{
    if (self->enemy)
        return;
    if (self->health <= 0)
        return;
    if (activator->flags & FL_NOTARGET)
        return;
    if (!(activator->client) && !(activator->monsterinfo.aiflags & AI_GOOD_GUY))
        return;

// delay reaction so if the monster is teleported, its sound is still heard
    self->enemy = activator;
    FoundTarget(self);
}


void monster_start_go(edict_t *self);


void monster_triggered_spawn(edict_t *self)
{
    self->s.origin[2] += 1;
    KillBox(self);

    self->solid = SOLID_BBOX;
    self->movetype = MOVETYPE_STEP;
    self->svflags &= ~SVF_NOCLIENT;
    self->air_finished_framenum = level.framenum + 12 * HZ;
    gi.linkentity(self);

    monster_start_go(self);

    if (self->enemy && !(self->spawnflags & 1) && !(self->enemy->flags & FL_NOTARGET)) {
        FoundTarget(self);
    } else {
        self->enemy = NULL;
    }
}

void monster_triggered_spawn_use(edict_t *self, edict_t *other, edict_t *activator)
{
    // we have a one frame delay here so we don't telefrag the guy who activated us
    self->think = monster_triggered_spawn;
    self->nextthink = level.framenum + 1;
    if (activator->client)
        self->enemy = activator;
    self->use = monster_use;
}

void monster_triggered_start(edict_t *self)
{
    self->solid = SOLID_NOT;
    self->movetype = MOVETYPE_NONE;
    self->svflags |= SVF_NOCLIENT;
    self->nextthink = 0;
    self->use = monster_triggered_spawn_use;
}


/*
================
monster_death_use

When a monster dies, it fires all of its targets with the current
enemy as activator.
================
*/
void monster_death_use(edict_t *self)
{
    self->flags &= ~(FL_FLY | FL_SWIM);
    self->monsterinfo.aiflags &= AI_GOOD_GUY;

    if (self->item) {
        Drop_Item(self, self->item);
        self->item = NULL;
    }

    // SPAQ
    while (self->monsterinfo.item_offset)
        M_DropItem(self);
    // SPAQ

    if (self->deathtarget)
        self->target = self->deathtarget;

    if (!self->target)
        return;

    G_UseTargets(self, self->enemy);
}

// SPAQ
static void SpawnChampionBodyguards(edict_t *self)
{
    int num_bodyguards = 2 + (rand() % (int)skill->value + 1);

    while (num_bodyguards)
    {
        for (int attempts = 0; attempts < 3; attempts++)
        {
            trace_t tr;
            float angle = random() * 360;

            vec3_t start, end, forward;
            VectorCopy(self->s.origin, start);

            AngleVectors((const vec_t[]) { 0, angle, 0 }, forward, NULL, NULL);

            VectorMA(start, (((self->mins[0] + self->mins[1]) - (self->maxs[0] + self->maxs[1])) / 2) + 24 + (random() * 120), forward, end);

            tr = gi.trace(start, self->mins, self->maxs, end, self, MASK_SHOT);

            if (tr.startsolid || tr.allsolid)
                continue;

            tr = gi.trace(tr.endpos, self->mins, self->maxs, tr.endpos, NULL, MASK_SHOT);

            if (tr.startsolid || tr.allsolid)
                continue;

            edict_t *bodyguard = G_Spawn();
            bodyguard->classname = self->classname;
            VectorCopy(tr.endpos, bodyguard->s.origin);
            bodyguard->s.angles[1] = angle;
            bodyguard->flags |= FL_NO_CHAMPS;
            ED_CallSpawn(bodyguard);

            break;
        }

        num_bodyguards--;
    }
}

static void SetupItems(edict_t *self)
{
    int ammo_count = (self->monsterinfo.champion ? 2 : 0) + rand() % 3;
    
    for (int i = 0; i < ammo_count; i++)
    {
        self->monsterinfo.items[self->monsterinfo.item_offset] = AMMO_FIRST + (rand() % (AMMO_COUNT + 2));

        // knives and grenades counted as ammo, so they drop more often
        if (self->monsterinfo.items[self->monsterinfo.item_offset] == AMMO_MAX)
            self->monsterinfo.items[self->monsterinfo.item_offset] = KNIFE_NUM;
        else if (self->monsterinfo.items[self->monsterinfo.item_offset] == AMMO_MAX + 1)
            self->monsterinfo.items[self->monsterinfo.item_offset] = GRENADE_NUM;

        self->monsterinfo.item_offset++;
    }

    qboolean spawn_weapon = self->monsterinfo.champion ? true : (random() < 0.1);

    if (spawn_weapon)
    {
        self->monsterinfo.items[self->monsterinfo.item_offset] = WEAPON_FIRST + (rand() % (WEAPON_COUNT - 2));

        // MK23 -> Dual
        if (self->monsterinfo.items[self->monsterinfo.item_offset] == MK23_NUM)
            self->monsterinfo.items[self->monsterinfo.item_offset] = DUAL_NUM;

        self->monsterinfo.item_offset++;
    }

    qboolean spawn_item = self->monsterinfo.champion ? true : (random() < 0.1);

    if (spawn_item)
        self->monsterinfo.items[self->monsterinfo.item_offset++] = ITEM_FIRST + (rand() % ITEM_COUNT);
}

static void CheckChampionSpawn(edict_t *self)
{
    static int offset = -1;

    if (!(self->flags & FL_NO_CHAMPS))
    {
        if (offset == -1)
            offset = rand() % (int) (6 - skill->value);

        if ((offset = ((offset + 1) % 4)) == 0)
        {
            self->monsterinfo.champion = rand() % CHAMPION_TOTAL;

            if (!self->solid)
                self->monsterinfo.champion &= ~CHAMPION_BODYGUARDS;
        }

        if (self->monsterinfo.champion)
        {
            if (MONSTER_IS_CHAMP(self, CHAMPION_STRONK))
                self->health *= (skill->value / 3) * 2;

            if (MONSTER_IS_CHAMP(self, CHAMPION_BODYGUARDS))
                SpawnChampionBodyguards(self);
        }
    }

    SetupItems(self);
}

#ifdef _WIN32
#include <malloc.h>
#endif

void CreateChampions(void)
{
#ifdef _WIN32
    int *champ_ids = _malloca(sizeof(int) * globals.num_edicts);
#else
    int champ_ids[globals.num_edicts];
#endif

    int num_monsters = 0;

    for (int i = game.maxclients + 1; i < globals.num_edicts; i++)
    {
        edict_t *e = &g_edicts[i];

        if (e->svflags & SVF_MONSTER)
            champ_ids[num_monsters++] = i;
    }

    for (int i = 0; i < num_monsters - 1; i++) 
    {
        size_t j = i + rand() / (RAND_MAX / (num_monsters - i) + 1);
        int t = champ_ids[j];
        champ_ids[j] = champ_ids[i];
        champ_ids[i] = t;
    }

    for (int i = 0; i < num_monsters; i++)
        CheckChampionSpawn(&g_edicts[champ_ids[i]]);
}

void M_DropItem(edict_t *self)
{
    if (!self->monsterinfo.item_offset)
        return;

    edict_t *dropped = Drop_Item(self, FindItemByNum(self->monsterinfo.items[--self->monsterinfo.item_offset]));

    vec3_t temp_angles = { 0, random() * 360, 0 };
    vec3_t forward;

	AngleVectors (temp_angles, forward, NULL, NULL);
	VectorScale (forward, 300, dropped->velocity);
    dropped->velocity[2] += 100;
}
// SPAQ

//============================================================================

qboolean monster_start(edict_t *self)
{
    if (deathmatch->value) {
        G_FreeEdict(self);
        return false;
    }

    if ((self->spawnflags & 4) && !(self->monsterinfo.aiflags & AI_GOOD_GUY)) {
        self->spawnflags &= ~4;
        self->spawnflags |= 1;
//      gi.dprintf("fixed spawnflags on %s at %s\n", self->classname, vtos(self->s.origin));
    }

    if (!(self->monsterinfo.aiflags & AI_GOOD_GUY))
        level.total_monsters++;

    self->nextthink = level.framenum + 1;
    self->svflags |= SVF_MONSTER;
    self->s.renderfx |= RF_FRAMELERP;
    self->takedamage = DAMAGE_AIM;
    self->air_finished_framenum = level.framenum + 12 * HZ;
    self->use = monster_use;
    self->max_health = self->health;
    self->clipmask = MASK_MONSTERSOLID;
    // SPAQ
    self->head_height = self->maxs[1];
    // SPAQ

    self->s.skinnum = 0;
    self->deadflag = DEAD_NO;
    self->svflags &= ~SVF_DEADMONSTER;

    if (!self->monsterinfo.checkattack)
        self->monsterinfo.checkattack = M_CheckAttack;
    VectorCopy(self->s.origin, self->s.old_origin);

    if (st.item) {
        self->item = FindItemByClassname(st.item);
        if (!self->item)
            gi.dprintf("%s at %s has bad item: %s\n", self->classname, vtos(self->s.origin), st.item);
    }

    // randomize what frame they start on
    if (self->monsterinfo.currentmove)
        self->s.frame = self->monsterinfo.currentmove->firstframe + (rand() % (self->monsterinfo.currentmove->lastframe - self->monsterinfo.currentmove->firstframe + 1));

    // SPAQ
    if (self->spawnflags & 2)
        self->flags |= FL_NO_CHAMPS;
    // SPAQ

    return true;
}

void monster_start_go(edict_t *self)
{
    vec3_t  v;

    if (self->health <= 0)
        return;

    // check for target to combat_point and change to combattarget
    if (self->target) {
        qboolean        notcombat;
        qboolean        fixup;
        edict_t     *target;

        target = NULL;
        notcombat = false;
        fixup = false;
        while ((target = G_Find(target, FOFS(targetname), self->target)) != NULL) {
            if (strcmp(target->classname, "point_combat") == 0) {
                self->combattarget = self->target;
                fixup = true;
            } else {
                notcombat = true;
            }
        }
        if (notcombat && self->combattarget)
            gi.dprintf("%s at %s has target with mixed types\n", self->classname, vtos(self->s.origin));
        if (fixup)
            self->target = NULL;
    }

    // validate combattarget
    if (self->combattarget) {
        edict_t     *target;

        target = NULL;
        while ((target = G_Find(target, FOFS(targetname), self->combattarget)) != NULL) {
            if (strcmp(target->classname, "point_combat") != 0) {
                gi.dprintf("%s at (%i %i %i) has a bad combattarget %s : %s at (%i %i %i)\n",
                           self->classname, (int)self->s.origin[0], (int)self->s.origin[1], (int)self->s.origin[2],
                           self->combattarget, target->classname, (int)target->s.origin[0], (int)target->s.origin[1],
                           (int)target->s.origin[2]);
            }
        }
    }

    if (self->target) {
        self->goalentity = self->movetarget = G_PickTarget(self->target);
        if (!self->movetarget) {
            gi.dprintf("%s can't find target %s at %s\n", self->classname, self->target, vtos(self->s.origin));
            self->target = NULL;
            self->monsterinfo.pause_framenum = INT_MAX;
            self->monsterinfo.stand(self);
        } else if (strcmp(self->movetarget->classname, "path_corner") == 0) {
            VectorSubtract(self->goalentity->s.origin, self->s.origin, v);
            self->ideal_yaw = self->s.angles[YAW] = vectoyaw(v);
            self->monsterinfo.walk(self);
            self->target = NULL;
        } else {
            self->goalentity = self->movetarget = NULL;
            self->monsterinfo.pause_framenum = INT_MAX;
            self->monsterinfo.stand(self);
        }
    } else {
        self->monsterinfo.pause_framenum = INT_MAX;
        self->monsterinfo.stand(self);
    }

    self->think = monster_think;
    self->nextthink = level.framenum + 1;
}


void walkmonster_start_go(edict_t *self)
{
    if (!(self->spawnflags & 2) && level.time < 1) {
        M_droptofloor(self);

        if (self->groundentity)
            if (!M_walkmove(self, 0, 0))
                gi.dprintf("%s in solid at %s\n", self->classname, vtos(self->s.origin));
    }

    if (!self->yaw_speed)
        self->yaw_speed = 20;
    self->viewheight = 25;

    monster_start_go(self);

    if (self->spawnflags & 2)
        monster_triggered_start(self);
}

void walkmonster_start(edict_t *self)
{
    self->think = walkmonster_start_go;
    monster_start(self);
}


void flymonster_start_go(edict_t *self)
{
    if (!M_walkmove(self, 0, 0))
        gi.dprintf("%s in solid at %s\n", self->classname, vtos(self->s.origin));

    if (!self->yaw_speed)
        self->yaw_speed = 10;
    self->viewheight = 25;

    monster_start_go(self);

    if (self->spawnflags & 2)
        monster_triggered_start(self);
}


void flymonster_start(edict_t *self)
{
    self->flags |= FL_FLY;
    self->think = flymonster_start_go;
    monster_start(self);
}


void swimmonster_start_go(edict_t *self)
{
    if (!self->yaw_speed)
        self->yaw_speed = 10;
    self->viewheight = 10;

    monster_start_go(self);

    if (self->spawnflags & 2)
        monster_triggered_start(self);
}

void swimmonster_start(edict_t *self)
{
    self->flags |= FL_SWIM;
    self->think = swimmonster_start_go;
    monster_start(self);
}
